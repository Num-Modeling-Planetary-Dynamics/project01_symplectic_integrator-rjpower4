from dataclasses import dataclass
import requests
import pathlib

OUTPUT_PATH = pathlib.Path(__file__).parent.resolve() / "kernels"


@dataclass
class Target:
    name: str
    url: str
    md5: str

    def output_path(self) -> pathlib.Path:
        """Return the path that the file should be saved to."""
        return OUTPUT_PATH / self.name

    def fetch(self, *, force=False):
        o_path = self.output_path()

        if force or not o_path.exists():
            with open(self.output_path(), "wb") as f:
                content = requests.get(self.url, stream=True).content
                f.write(content)
        elif o_path.exists():
            print(f"WARNING: Skipping {o_path} as it exists (set force=True to overwrite)")


TARGETS = [
    Target(
        "naif0012.tls",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls",
        "25a2fff30b0dedb4d76c06727b1895b1",
    ),
    Target(
        "gm_de431.tpc",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc",
        "6445f12003d1effcb432ea2671a51f18",
    ),
    Target(
        "de440.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
        "c9d581bfd84209dbeee8b1583939b148",
    ),
    Target(
        "nep097.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/nep097.bsp",
        "937e05462aac5e34cb5cfddae6c99199",
    ),
    Target(
        "plu058.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/plu058.bsp",
        "f1173cd80bd13feb91d0e29487292bbf",
    ),
]

for target in TARGETS:
    if not OUTPUT_PATH.is_dir():
        OUTPUT_PATH.mkdir()
    target.fetch()
