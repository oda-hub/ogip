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

    S = ogip.core.open_something(fn)

    print(S.to_long_string())


@spec.command()
@click.argument("input_fn")
@click.argument("output_fn")
def fix(input_fn, output_fn):
    """read and write, ensuring the format"""

    S = ogip.spec.Spectrum.from_file_name(input_fn)

    print(S.to_long_string())

    S.to_fits(output_fn)

if __name__ == "__main__":
    cli()
