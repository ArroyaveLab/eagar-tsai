"""Setup script for the eagar-tsai C extension."""

import sys

from setuptools import Extension, setup

extra_args: list[str]
if sys.platform == "win32":
    extra_args = ["/O2"]
else:
    extra_args = ["-O3", "-ffast-math"]

setup(
    ext_modules=[
        Extension(
            "eagar_tsai._integrand_ext",
            sources=["src/eagar_tsai/_integrand_ext.c"],
            extra_compile_args=extra_args,
        )
    ]
)
