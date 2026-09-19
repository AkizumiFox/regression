r"""Numbers computed in code, recorded so tools/check_numbers.py can confirm the book still quotes them.

    gen = Generated("ch06", "fwl")              # prefix defaults to the script name
    gen.num("beta", 0.4123, digits=3)           # key ch06:fwl:beta
    gen.write()
"""
from .paths import generated_path


class Generated:
    def __init__(self, chapter: str, name: str, prefix: str | None = None):
        self.chapter, self.name = chapter, name
        self.prefix = prefix or name
        self.items: dict[str, str] = {}

    def text(self, key: str, value: str):
        self.items[key] = value
        return value

    def num(self, key: str, value: float, digits: int = 3, sci: bool = False):
        if sci:
            mant, exp = f"{value:.{digits}e}".split("e")
            s = f"{mant}\\times 10^{{{int(exp)}}}"
            s = f"\\ensuremath{{{s}}}"
        else:
            s = f"{value:.{digits}f}"
            if s.startswith("-"):
                s = "\\ensuremath{-}" + s[1:]
        return self.text(key, s)

    def int(self, key: str, value: int):
        return self.text(key, f"{int(value)}")

    def write(self):
        import json
        path = generated_path(self.chapter, self.name)
        values = {f"{self.chapter}:{self.prefix}:{k}": v for k, v in self.items.items()}
        path.write_text(json.dumps(values, indent=1, ensure_ascii=False) + "\n")
        print(f"wrote {path.relative_to(path.parents[3])} ({len(self.items)} values)")
