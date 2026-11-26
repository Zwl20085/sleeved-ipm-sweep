"""Parser for the plain-text AEDT project format (``$begin``/``$end`` blocks).

The file is a tree of named blocks. Every non-block line is a field of one of
three shapes::

    key=value                     (assign)
    Name(arg, arg, key=value)     (call)
    Name[n: item, item]           (array)

Values are parsed lazily with :func:`parse_value`; the raw text is always kept
so nothing is lost.  Only the ``AnsoftProject`` text part is parsed; the binary
``AllReferencedFilesForProject`` tail is located separately (see
:func:`read_native_files`).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Iterator, Optional, Union


class AedtParseError(ValueError):
    """Raised when a line or value cannot be parsed."""


@dataclass
class Field:
    key: str
    kind: str  # 'assign' | 'call' | 'array' | 'bare'
    raw: str  # raw value text (for call/array: the inside of the brackets)
    line: int

    def value(self) -> Any:
        if self.kind == "assign":
            return parse_value(self.raw)
        if self.kind == "call":
            return parse_args(self.raw)
        if self.kind == "array":
            return parse_array_body(self.raw)
        return self.raw


@dataclass
class Block:
    name: str
    line: int
    end_line: int = -1
    items: list[Union["Block", Field]] = field(default_factory=list)

    # ---- navigation -------------------------------------------------
    def blocks(self, name: Optional[str] = None) -> list["Block"]:
        return [b for b in self.items if isinstance(b, Block) and (name is None or b.name == name)]

    def block(self, name: str) -> "Block":
        found = self.blocks(name)
        if len(found) != 1:
            raise KeyError(f"block {name!r} occurs {len(found)} times in {self.name!r} (line {self.line})")
        return found[0]

    def opt_block(self, name: str) -> Optional["Block"]:
        found = self.blocks(name)
        if len(found) > 1:
            raise KeyError(f"block {name!r} occurs {len(found)} times in {self.name!r} (line {self.line})")
        return found[0] if found else None

    def fields(self, key: Optional[str] = None) -> list[Field]:
        return [f for f in self.items if isinstance(f, Field) and (key is None or f.key == key)]

    def field(self, key: str) -> Field:
        found = self.fields(key)
        if len(found) != 1:
            raise KeyError(f"field {key!r} occurs {len(found)} times in {self.name!r} (line {self.line})")
        return found[0]

    def get(self, key: str, default: Any = None) -> Any:
        found = self.fields(key)
        if not found:
            return default
        return found[-1].value()

    def __getitem__(self, key: str) -> Any:
        return self.field(key).value()

    def has(self, key: str) -> bool:
        return bool(self.fields(key))

    def walk(self) -> Iterator["Block"]:
        yield self
        for b in self.blocks():
            yield from b.walk()

    def to_dict(self) -> dict:
        """Plain nested dict (repeated keys become lists)."""
        out: dict[str, Any] = {}
        for it in self.items:
            if isinstance(it, Block):
                key, val = it.name, it.to_dict()
            else:
                key = it.key
                try:
                    val = it.value()
                except AedtParseError:
                    val = it.raw
            if key in out:
                if not isinstance(out[key], list) or not getattr(out[key], "_multi", False):
                    lst = _MultiList([out[key]])
                    out[key] = lst
                out[key].append(val)
            else:
                out[key] = val
        return {k: (list(v) if isinstance(v, _MultiList) else v) for k, v in out.items()}


class _MultiList(list):
    _multi = True


# ---------------------------------------------------------------------------
# line level
_KEY_RE = re.compile(r"^(?:'((?:[^'\\]|\\.)*)'|([A-Za-z_0-9][^=(\[\s]*))")


def parse_line(text: str, lineno: int) -> Field:
    m = _KEY_RE.match(text)
    if not m:
        return Field(key="", kind="bare", raw=text, line=lineno)
    key = m.group(1) if m.group(1) is not None else m.group(2)
    rest = text[m.end():]
    if rest.startswith("="):
        return Field(key=key, kind="assign", raw=rest[1:], line=lineno)
    if rest.startswith("(") and rest.endswith(")"):
        return Field(key=key, kind="call", raw=rest[1:-1], line=lineno)
    if rest.startswith("[") and rest.endswith("]"):
        return Field(key=key, kind="array", raw=rest[1:-1], line=lineno)
    return Field(key=key, kind="bare", raw=text, line=lineno)


def parse_lines(lines: list[str], start: int, stop: int, first_lineno: int = 1) -> Block:
    """Parse ``lines[start:stop]`` which must contain exactly one top block."""
    root = Block(name="<root>", line=first_lineno + start)
    stack = [root]
    for i in range(start, stop):
        text = lines[i].strip()
        lineno = first_lineno + i
        if not text:
            continue
        if text.startswith("$begin '"):
            b = Block(name=text[8:-1], line=lineno)
            stack[-1].items.append(b)
            stack.append(b)
        elif text.startswith("$end '"):
            name = text[6:-1]
            if len(stack) < 2 or stack[-1].name != name:
                raise AedtParseError(f"line {lineno}: unbalanced $end {name!r} (open: {stack[-1].name!r})")
            stack[-1].end_line = lineno
            stack.pop()
        else:
            stack[-1].items.append(parse_line(text, lineno))
    if len(stack) != 1:
        raise AedtParseError(f"unclosed block {stack[-1].name!r} opened at line {stack[-1].line}")
    return root


# ---------------------------------------------------------------------------
# value level
_TOKEN_RE = re.compile(
    r"\s*(?:(?P<str>'(?:[^'\\]|\\.)*')|(?P<num>[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?(?![A-Za-z_]))"
    r"|(?P<punct>[(),=\[\]:])|(?P<word>[^\s(),=\[\]:']+))"
)


def _tokenize(text: str) -> list[tuple[str, str]]:
    pos, out = 0, []
    while pos < len(text):
        m = _TOKEN_RE.match(text, pos)
        if not m or m.end() == pos:
            if text[pos:].strip() == "":
                break
            raise AedtParseError(f"cannot tokenize value at {text[pos:pos+40]!r} in {text[:80]!r}")
        for kind in ("str", "num", "punct", "word"):
            if m.group(kind) is not None:
                out.append((kind, m.group(kind)))
                break
        pos = m.end()
    return out


@dataclass
class Call:
    name: str
    args: list[Any]
    kwargs: dict[str, Any]


def _atom(kind: str, tok: str) -> Any:
    if kind == "str":
        return tok[1:-1].replace("\\'", "'")
    if kind == "num":
        return float(tok) if any(c in tok for c in ".eE") else int(tok)
    if tok == "true":
        return True
    if tok == "false":
        return False
    return tok


class _ArgParser:
    def __init__(self, text: str):
        self.text = text
        self.toks = _tokenize(text)
        self.i = 0

    def peek(self) -> Optional[tuple[str, str]]:
        return self.toks[self.i] if self.i < len(self.toks) else None

    def take(self) -> tuple[str, str]:
        if self.i >= len(self.toks):
            raise AedtParseError(f"unexpected end of value {self.text[:80]!r}")
        t = self.toks[self.i]
        self.i += 1
        return t

    def value(self) -> Any:
        kind, tok = self.take()
        nxt = self.peek()
        if kind == "punct" and tok == "(":
            items = self.items(")")
            return [v for v in items]
        if kind == "punct" and tok == "[":
            return self.array_rest()
        if kind == "punct":
            raise AedtParseError(f"unexpected {tok!r} in {self.text[:80]!r}")
        if nxt == ("punct", "(") and kind in ("word", "str"):
            self.take()
            args, kwargs = self.args_until(")")
            return Call(name=tok if kind == "word" else tok[1:-1], args=args, kwargs=kwargs)
        if nxt == ("punct", "[") and kind == "word":
            self.take()
            return Call(name=tok, args=self.array_rest(), kwargs={})
        return _atom(kind, tok)

    def array_rest(self) -> list[Any]:
        # '[n: a, b]' ; the leading count is checked
        kind, tok = self.take()
        if kind != "num":
            raise AedtParseError(f"array without count in {self.text[:80]!r}")
        n = int(tok)
        if self.take() != ("punct", ":"):
            raise AedtParseError(f"array without ':' in {self.text[:80]!r}")
        items = self.items("]")
        if len(items) != n:
            raise AedtParseError(f"array count {n} != {len(items)} items in {self.text[:80]!r}")
        return items

    def items(self, close: str) -> list[Any]:
        args, kwargs = self.args_until(close)
        if kwargs:
            return args + [kwargs]
        return args

    def args_until(self, close: Optional[str]) -> tuple[list[Any], dict[str, Any]]:
        args: list[Any] = []
        kwargs: dict[str, Any] = {}
        while True:
            t = self.peek()
            if t is None:
                if close is None:
                    return args, kwargs
                raise AedtParseError(f"missing {close!r} in {self.text[:80]!r}")
            if t == ("punct", close):
                self.take()
                return args, kwargs
            if t == ("punct", ","):
                self.take()
                continue
            # key=value ?
            if self.i + 1 < len(self.toks) and self.toks[self.i + 1] == ("punct", "=") and t[0] in ("word", "str"):
                key = _atom(*t)
                self.i += 2
                kwargs[str(key)] = self.value()
            else:
                args.append(self.value())


def parse_value(raw: str) -> Any:
    p = _ArgParser(raw)
    v = p.value()
    if p.peek() is not None:
        # e.g. a bare word list -- keep raw
        raise AedtParseError(f"trailing tokens in value {raw[:80]!r}")
    return v


def parse_args(raw: str) -> list[Any]:
    """Arguments of a call line ``Name(...)``; kwargs appended as a trailing dict."""
    p = _ArgParser(raw)
    args, kwargs = p.args_until(None)
    return args + ([kwargs] if kwargs else [])


def parse_array_body(raw: str) -> list[Any]:
    p = _ArgParser("[" + raw + "]")
    p.take()
    return p.array_rest()


# ---------------------------------------------------------------------------
# file level
@dataclass
class AedtFile:
    path: str
    lines: list[str]
    text_end: int  # index (0-based, exclusive) of the end of 'AnsoftProject'

    @classmethod
    def load(cls, path: str) -> "AedtFile":
        with open(path, "rb") as fh:
            data = fh.read()
        marker = data.find(b"\n$begin 'AllReferencedFilesForProject'")
        text = data[: marker + 1] if marker >= 0 else data
        lines = text.decode("utf-8", errors="strict").split("\n")
        # lines[] is 0-based; line number = index + 1
        end = None
        for i in range(len(lines) - 1, -1, -1):
            if lines[i].strip() == "$end 'AnsoftProject'":
                end = i + 1
                break
        if end is None:
            raise AedtParseError(f"{path}: no $end 'AnsoftProject'")
        return cls(path=path, lines=lines, text_end=end)

    def find_block_end(self, start: int) -> int:
        """Given 0-based index of a ``$begin`` line, return index after its ``$end``."""
        head = self.lines[start]
        indent = head[: len(head) - len(head.lstrip("\t"))]
        name = head.strip()[8:-1]
        target = f"{indent}$end '{name}'"
        for i in range(start + 1, self.text_end):
            if self.lines[i].rstrip("\r") == target:
                return i + 1
        raise AedtParseError(f"{self.path}: no end for block {name!r} at line {start+1}")

    def parse_block_at(self, start: int) -> Block:
        stop = self.find_block_end(start)
        root = parse_lines(self.lines, start, stop, first_lineno=1)
        return root.items[0]  # type: ignore[return-value]

    def top_level_starts(self, depth: int, name: str) -> list[int]:
        prefix = "\t" * depth + f"$begin '{name}'"
        return [i for i in range(self.text_end) if self.lines[i].rstrip("\r") == prefix]

    def design_blocks(self) -> dict[str, Block]:
        out: dict[str, Block] = {}
        for kind in ("Maxwell2DModel", "Maxwell3DModel"):
            for s in self.top_level_starts(1, kind):
                b = self.parse_block_at(s)
                out[b["Name"]] = b
        return out

    def project_block(self, name: str) -> Block:
        starts = self.top_level_starts(1, name)
        if len(starts) != 1:
            raise AedtParseError(f"{self.path}: block {name!r} found {len(starts)} times at depth 1")
        return self.parse_block_at(starts[0])


def read_native_files(path: str) -> dict[str, bytes]:
    """Embedded native geometry files (Parasolid ``.x_b``) by stored path."""
    with open(path, "rb") as fh:
        data = fh.read()
    out: dict[str, bytes] = {}
    pat = re.compile(rb"\$begin 'x_b'\r?\n([^\r\n]+)\r?\nBIN(\d{12})")
    for m in pat.finditer(data):
        name = m.group(1).decode()
        n = int(m.group(2))
        start = m.end()
        out[name] = data[start: start + n]
    return out
