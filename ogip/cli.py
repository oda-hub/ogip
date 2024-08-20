import click

import ogip.spec
import ogip.core


@click.group()
def cli():
    pass


@cli.group()
def spec():
    pass


@spec.command()
@click.argument("FN")
def inspect(fn):
    s = ogip.core.open_something(fn)

    print(s.to_long_string())


@spec.command()
@click.argument("input_fn")
@click.argument("output_fn")
def fix(input_fn, output_fn):
    """read and write, ensuring the format"""

    s = ogip.spec.Spectrum.from_file_name(input_fn)

    print(s.to_long_string())

    s.to_fits(output_fn)


if __name__ == "__main__":
    cli()
