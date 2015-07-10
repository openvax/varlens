import matplotlib
from matplotlib import pyplot
import os

from matplotlib.backends.backend_pdf import PdfPages

def configure_matplotlib(force_agg=False, interactive=False):
    if not os.environ.get('DISPLAY') or force_agg:
        matplotlib.rcParams['backend'] = 'Agg'
    matplotlib.interactive(interactive)

def write_figures(filename, figures):
    out_extension = filename.split('.')[-1]
    if out_extension == 'pdf':
        if not figures:
            raise ValueError("No figures to write.")

        pages = PdfPages(filename)
        for figure in figures:
            pages.savefig(figure)
        pages.close()
        print("Wrote: %s." % filename)
            
    elif out_extension in ('png', 'eps', 'ps', 'svg'):
        if not figures:
            raise ValueError("No figures to write.")

        (fig,) = figures
        fig.savefig(filename)
        print("Wrote: %s." % filename)

    elif filename == 'display':
        pyplot.show()
    elif filename == 'none':
        pass
    else:
        print("Unsupported output: %s" % filename)
