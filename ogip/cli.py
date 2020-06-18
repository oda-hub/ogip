import click

import ogip.spec

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
@click.argument("input_fn", help="read and write, ensuring the format")
@click.argument("output_fn", help="read and write, ensuring the format")
def fix(input_fn, output_fn):

    S = ogip.spec.Spectrum.from_file_name(input_fn)

    print(S.to_long_string())

    S.to_fits(output_fn)

if __name__ == "__main__":
    cli()
