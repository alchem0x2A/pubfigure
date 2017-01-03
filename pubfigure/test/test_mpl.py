from pubfigure.FigureCollection import *
from scipy import sin, cos, tan, pi
import numpy

def cot(x):
    return 1/tan(x)

def plot_fun(fig, i):
    funcs=[sin, cos, tan, cot]
    x = numpy.linspace(-2*pi, 2*pi, 100)
    y = funcs[i](x)
    ax = fig.add_subplot(111)
    ax.plot(x,y)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(funcs[i].__name__+r"($x$)")
    if i<2:
        ax.set_ylim(-1.0, 1.0)
    else:
        ax.set_ylim(-10, 10)
    ax.set_xticks(numpy.array([-2, -1.5, -1, -0.5, 0,
                               0.5, 1, 1.5, 2])*pi)
    ax.set_xticklabels([r"$-2\pi$", r"$-3\pi/2$",
                        r"$-\pi$", r"$-\pi/2$",
                        "0",
                        r"$\pi/2$", r"$\pi$",
                        r"$3\pi/2$", r"$2\pi$"])
    fig.tight_layout(pad=0)

fc = FigureCollection(figure_style="nature",
                      collection_style="nature",
                      pagesize=(6.0, 4.5),
                      col=2,
                      row=2)

for i in range(4):
    fig, num = fc.add_figure()
    fig.set_plot_func(plot_fun, i)
    # print (fig.x_, fig.y_, fig.w_, fig.h_)
    # print(fig.fig_param["annotation.location"])
    print(fig.get_x_text(), fig.get_y_text())
    print (fig.get_w_fig(), fig.get_h_fig())
    # print(fig.fig_param["figure.lpad"], fig.fig_param["figure.rpad"],
          # fig.fig_param["figure.tpad"], fig.fig_param["figure.bpad"])

fc.save_all("./test_mpl/test_mpl.pdf", outline=True)
