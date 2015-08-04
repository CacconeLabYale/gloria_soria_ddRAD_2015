
import click



# this is the command function

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('--test-opt',
    show_default=True,
    default='test-it',
    help="testing tests the test tests")
@click.argument('test_arg', help='arg help')
def test(test_opt):
    pass


# @click.option('--distance-bin', default=50, show_default=True, help="How wide do you want the bin window?")
# @click.argument('ld_path',type=click.Path(exists=True),help="Path to the table file created by the LD calculation program.")
# @click.argument('out_path', type=click.Path(), help="Path to where you want to save the results pickle.")
# @click.option('--count', default=1, help='number of greetings')