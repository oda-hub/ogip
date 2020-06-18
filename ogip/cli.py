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

    S = ogip.spec.Spectrum.from_file_name(fn)

    print(S.to_long_string())

if __name__ == "__main__":
    cli()
