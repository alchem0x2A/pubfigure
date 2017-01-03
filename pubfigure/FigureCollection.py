"""
:class: `FigureCollector`
    Collect the subfigures generated from matplotlib.
    The idea is to align the subfigures and unify the style,
    especially the fontsize. 
    `Tikz` package is required to generate the final page.
"""
import matplotlib
import matplotlib.style
from matplotlib.cbook import Bunch
import pkg_resources
import warnings
from matplotlib.figure import Figure
from matplotlib.backends import backend_agg
import shutil
import os
import os.path
import jinja2
import subprocess


package_name = __name__.split(".")[0]
_append_file_style_folder = lambda s, end: "/".join(("styles", s+end))
style_dir_name = pkg_resources.resource_filename(package_name, "styles")
template_dir_name = pkg_resources.resource_filename(package_name, "templates")


class TemplateFiller(object):
    def __init__(self,
                 temp_dir=template_dir_name,
                 template="fc_tikz.tex"):
        """
        Use filler to fill in template
        TODO: other engines?
        """
        self.jinja_env = jinja2.Environment(
            block_start_string='\BLOCK{',
            block_end_string='}',
            variable_start_string='\VAR{',
            variable_end_string='}',
            comment_start_string='\#{',
            comment_end_string='}',
            line_statement_prefix='%%',
            line_comment_prefix='%#',
            trim_blocks=True,
            autoescape=False,
            loader=jinja2.FileSystemLoader(os.path.abspath(temp_dir))
        )
        self.template = self.jinja_env.get_template(template)

    def render(self, filename=None, **argv):
        s = self.template.render(**argv)
        if filename is None:
            return s            # Only returns if no file saves
        else:
            try:
                with open(filename, "w") as f:
                    f.write(s)
                return filename
            except PermissionError:
                raise


# Legal parameter keys for FigureCollection
fc_param_keys = ("page.width",
                 "page.height",
                 "engine",
                 "figure.lpad",
                 "figure.rpad",
                 "figure.tpad",
                 "figure.bpad",
                 "annotation.font",
                 "annotation.style",
                 "annotation.size",
                 "annotation.weight",
                 "annotation.fc",
                 "annotation.bc",
                 "annotation.alpha",
                 "annotation.location")

default_fc_param = {"page.width": 8.0,
                    "page.height": 6.0,
                    "engine": "Tikz",
                    "figure.lpad": 0.1,
                    "figure.rpad": 0.1,
                    "figure.tpad": 0.1,
                    "figure.bpad": 0.1,
                    "annotation.font": "Arial",
                    "annotation.style": None,
                    "annotation.size": None,
                    "annotation.weight": None,
                    "annotation.fc": None,
                    "annotation.bc": None,
                    "annotation.alpha": None,
                    "annotation.location": None}


def _set_param_value(val):
    """
    The initial value passed is a `str`.
    Use "," to separate values and convert digits to float
    """
    def _to_digit(s):
        if s.strip() == "None":
            return None
        try:
            int(s)
            return int(s)
        except ValueError:
            try:
                float(s)
                return float(s)
            except ValueError:
                return s

    unpack = val.strip().split(",")
    if len(unpack) == 1:
        return _to_digit(unpack[0])
    else:
        return tuple(_to_digit(_v) for _v in unpack)


def _fc_params_in_file(fname, fail_on_error=False):
    """
    Like the function `matplotlib._rc_params_in_file`,
    but handles parameters for FigureCollection.
    """
    cnt = 0
    rc_temp = {}
    with open(fname, "r") as fd:
        try:
            for line in fd:
                cnt += 1
                strippedline = line.split('#', 1)[0].strip()
                if not strippedline:
                    continue
                tup = strippedline.split(':', 1)
                if len(tup) != 2:
                    error_details = "%d, $d, $s" % (cnt, line, fname)
                    warnings.warn('Illegal %s' % error_details)
                    continue
                key, val = tup
                key = key.strip()
                val = val.strip()
                if key in rc_temp:
                    warnings.warn('Duplicate key in file "%s", line #%d' %
                                  (fname, cnt))
                rc_temp[key] = (val, line, cnt)
        except UnicodeDecodeError:
            warnings.warn(
                ('Cannot decode configuration file for FigureCollection %s'
                 'check LANG and LC_* variables')
                % fname)
            raise

    config = {}

    for key in fc_param_keys:
        if key in rc_temp:
            val, line, cnt = rc_temp.pop(key)
            if fail_on_error:
                config[key] = _set_param_value(val)  # try to convert to proper type or raise
            else:
                try:
                    config[key] = _set_param_value(val)  # try to convert to proper type or skip
                except Exception as msg:
                    # error_details = _error_details_fmt % (cnt, line, fname)
                    warnings.warn('Bad val "%s" \n%s' %
                                  (val, msg))
        else:                   # Use default parameter
            config[key] = default_fc_param[key]

    return config


def _gen_combination(left, right, case="lower"):
    from string import ascii_lowercase, ascii_uppercase
    from itertools import product
    if case == "num":
        i = 1
        while True:
            yield "".join((left, str(i), right))
            i += 1
    elif case == "lower":
        size = 1
        while True:
            for s in product(ascii_lowercase, repeat=size):
                yield "".join((left, ) + s + (right, ))
            size += 1
    elif case == "upper":
        size = 1
        while True:
            for s in product(ascii_uppercase, repeat=size):
                yield "".join((left, ) + s + (right, ))
            size += 1
    else:
        raise ValueError("The case can only be num, lower or upper!")


class FigureCollection(object):
    def _load_style(self, default_style):
        # If the default style is user-defined or contained in the mpl
        # library, use the setting for default figure plot.
        if default_style is None:
            return

        if default_style in matplotlib.style.available:
            matplotlib.style.use(default_style)
        elif default_style.endswith(".mplstyle"):  # The user may use a defined mplstyle sheet, with absolute path!
            try:
                matplotlib.style.use(default_style) # TODO: add assert for non-existing file
            except FileNotFoundError:
                warnings.warn(("Warning: Unable to find the requested mplfigure style [%s],"
                              " will fallback to original rcParams. \n")
                              % default_style)
        elif pkg_resources.resource_exists(package_name,
                                           _append_file_style_folder(default_style, ".mplstyle")):
            # Use the style-sheet with this package
            matplotlib.style.use(pkg_resources.resource_filename(package_name,
                                                                 _append_file_style_folder(default_style, ".mplstyle")))
        else:
            warnings.warn("Warning: Unable to find the requested figure style [%s],"\
                          " will fallback to original rcParams. \n" % default_style)

    def _load_collection_style(self, collection_style):
        """
        TODO: add support for user-based collection style file!
        """
        self.fc_param = {}
        if collection_style == None:
            self.fc_param = default_fc_param
            return
            
        if collection_style.endswith(".fc"):  # Assume the user provides a valid path to name
            try:
                self.fc_param = _fc_params_in_file(collection_style)
            except FileNotFoundError:
                warnings.warn("The file %s cannot be found, Use default params for FigureCollection" % collection_style)
                self.fc_param = default_fc_param
        elif pkg_resources.resource_exists(package_name,
                                           _append_file_style_folder(collection_style, ".fc")):
            # Use the style-sheet with this package
            self.fc_param = _fc_params_in_file(pkg_resources.resource_filename(package_name,
                                                _append_file_style_folder(collection_style, ".fc")))
        else:
            warnings.warn(("Warning: Unable to find the requested FigureCollection style [%s],"
                          " will fallback to original rcParams. \n")
                          % collection_style)                 # The file is within the user folder

    def _treat_annotation_style(self, ann_str):
        if ann_str is None:
            return None
        else:
            _ann_dic = ann_str.strip().split("$")
            if len(_ann_dic) is not 3:
                raise IndexError("The style of annotation should be chained by $$")
            elif _ann_dic[1] == "a":
                gen = _gen_combination(_ann_dic[0], _ann_dic[2], "lower")
                return gen
            elif _ann_dic[1] == "A":
                gen = _gen_combination(_ann_dic[0], _ann_dic[2], "upper")
                return gen
            elif _ann_dic[1] == "1":
                gen = _gen_combination(_ann_dic[0], _ann_dic[2], "num")
                return gen
            else:
                raise ValueError("Only accept num, a or A as the style character!")

    def __init__(self,
                 pagesize=None,
                 figure_style=None,
                 collection_style="classic",
                 engine="Tikz",
                 annotation_style=None,
                 row=None,
                 col=None):
        """
        `pagesize` gives the size in ABSOLUTE inches
        `row` and `col` defines the rows and column numbers
        All class members with length unit, `self.w`, `self.h`
        are in absolute units;
        All parameters in fc_params are relative units.
        """
        self._mpl_rc_copy = matplotlib.rcParams.copy()
        # Load style sheet
        self._load_style(figure_style)
        self._load_collection_style(collection_style)
        # Overwrite the fc settings by user setting
        if pagesize is not None:
            try:
                self.fc_param["page.width"] = float(pagesize[0])
                self.fc_param["page.height"] = float(pagesize[1])
            except (ValueError, IndexError, TypeError):
                warnings.warn("pagesize must take a tuple of floats! \n")
                raise

        self.w = self.fc_param["page.width"]
        self.h = self.fc_param["page.height"]

        if annotation_style is not None:
            self.fc_param["annotation.style"] = annotation_style

        # Generator of annotation string. Use `next(self.generator)`
        # to generate a new style
        if self.fc_param["annotation.style"] is None:
            self.generator = None
        else:
            self.generator = self._treat_annotation_style(self.fc_param["annotation.style"])

        if engine is not None:
            self.fc_param["engine"] = engine

        self.use_grid = False
        if row is not None and col is not None:
            self.add_grid(col, row)

        self.figure_list = []
        self.cnt = 0            # count for total number of figures

    def add_grid(self, col, row):
        """
        Use grid to define the spacing of subfigures
        """
        self.use_grid = True
        if col > 0 and row > 0:
            self.col = col
            self.row = row
            self.sp_col = self.w / self.col
            self.sp_row = self.h / self.row
        else:
            raise ValueError("col and row must be positive!")

    def _treat_grid(self, loc):
        """
        Takes care of the location if `self.use_grid` is True
        The indices should start from 0
        """
        if type(loc) is int:       # loc takes exactly the i-th grid, row-order
            if loc >= self.col*self.row:
                raise ValueError("Current number is larger than total grids!")
            x_ = self.sp_col*(loc % self.col)
            y_ = self.sp_row*(loc // self.col)
            return (x_, y_, self.sp_col, self.sp_row)
        else:
            if len(loc) is 2:     # Use only 1 grid
                for d in loc:
                    if (int(d) is not d) or (d < 0):
                        raise ValueError("Must take positive intergers as loc for grid! when loc has 2 members")
                if (loc[0] >= self.col) or (loc[1] >= self.row):
                    raise ValueError("Current number is larger than total grids!")
                x_ = self.sp_col*loc[0]
                y_ = self.sp_row*loc[1]
                return (x_, y_, self.sp_col, self.sp_row)
            elif len(loc) is 4:  # Grid with grid span
                if False not in [0.0 <= i <= 1.0 for i in loc]:  # absolute loc value is used!
                    if (loc[0]+loc[2] > 1.0) or (loc[1]+loc[3] > 1.0):
                        warnings.warn("The figure use is larger than the whole page!")
                    x_ = self.w*loc[0]
                    y_ = self.h*loc[1]
                    w_ = self.w*loc[2]
                    h_ = self.h*loc[3]
                    return (x_, y_, w_, h_)
                else:
                    for d in loc:
                        if (int(d) is not d) or (d < 0):
                            raise ValueError("Must take positive intergers as loc for grid!")
                if (loc[0] >= self.col) or (loc[1] >= self.row):
                    raise ValueError("Current number is larger than total grids!")
                x_ = self.sp_col*loc[0]
                y_ = self.sp_row*loc[1]
                w_ = self.sp_col*loc[2]
                h_ = self.sp_row*loc[3]
                return (x_, y_, w_, h_)
            else:               # Too many members to pack!
                raise ValueError("In valid number of loc parameters!")

    def add_figure(self, loc=None, label=True, fig_file=None, **argv):
        # TODO add grid support
        """
        `label` determines whether to add the annotation label or not.
        `**argv` the user-defined fig_param
        User can also pass the `fig` instance to the subfigure
        """
        if self.use_grid is False:
            if loc is None:
                raise ValueError("Must define grids first!")
            if (type(loc) is not tuple) and (type(loc) is not list):
                raise TypeError("The absolute loc must be iterable!")
            if False in [0.0 <= i <= 1.0 for i in loc]:  # The values of loc should be within 1 and 0!!
                raise ValueError("The absolute loc value should be within 0 and 1!")
            loc_ = loc          # If grid is not used, the user must give a absolute position
        else:
            if loc is None:
                loc = self.cnt  # assume that the grids are added sequentially
            loc_ = self._treat_grid(loc)
        subfig_ = SubFigure(loc_, fig_file=fig_file,
                            parent_param=self.fc_param, **argv)
        if label is True:
            if self.generator is not None:
                s = next(self.generator)
                subfig_.set_annotation_text(s)
            else:
                s = None
        else:
            s = None
        b_ = Bunch(subfig=subfig_,
                   cnt=self.cnt,
                   label=label,
                   text=s,
                   file_name=None)
        self.figure_list.append(b_)
        self.cnt += 1           # In case the adding for figure failed, esp. in commandline mode
        return subfig_, b_.cnt       # User can use the returned instance to call the plot function

    def get_figure(self, num):
        if num > len(self.figure_list):
            raise ValueError("The requested figure number %d is larger than total figure numbers!" % num)
        return self.figure_list[num].subfig

    def save_one(self, num, filename):
        fig = self.get_figure(num)
        fig.save_fig(filename)

    def save_all(self, filename, save_figs=True,
                 export_tikz=True, outline=False):
        """
        Save the FigureCollection as a pdf file.
        `save_figs` controls if all subfigures are saved
        `export_tikz` exports a latex file with tikz
        TODO: save as other forms, e.g. png
        TODO: use other engines other than tikz
        TODO: assymetrically export figures
        """
        dir_name = os.path.dirname(os.path.abspath(filename))
        f_name = os.path.basename(os.path.abspath(filename))
        # TODO: add other file type support!
        if f_name.endswith(".pdf") is False:
            warnings.warn("Only pdf is supported!")
        # f_base = f_name.strip().split(".")[0]  # Sorry but dot it better not used in the file
        f_base = os.path.splitext(f_name)[0]
        dir_tmp_name = os.path.join(dir_name, f_base+"_tmp/")
        if os.path.exists(dir_tmp_name):
            shutil.rmtree(dir_tmp_name)  # I suppose you do not need the temporary files
        while os.path.exists(dir_tmp_name):  # Seems the rmtree takes sometime
            pass
        os.makedirs(dir_tmp_name)
        while os.path.exists(dir_tmp_name) is False:
            pass
        for i in range(len(self.figure_list)):
            fig_ = self.get_figure(i)
            if fig_.fig_type is "mpl":
                tmp_name = os.path.join(dir_tmp_name, "%d.pdf" % i)
            else:
                extension = os.path.splitext(fig_.fig)[-1]  # Nicer split of file extension
                tmp_name = os.path.join(dir_tmp_name, format(i)+extension)
            print(tmp_name)
            fig_.save_fig(tmp_name)
            self.figure_list[i].file_name = tmp_name

        # TODO: export using tikz here
        self._export(dir_name, f_base, outline)

        if save_figs is False:
            shutil.rmtree(dir_tmp_name)

        if export_tikz is True:
            # TODO: export tikz here
            pass

    def _export(self, dir_name, f_base, outline):
        """
        Export figures using the export engine.
        Output: .tex and .pdf (for tikz engine)
        """
        if self.fc_param["engine"] is "Tikz":
            file_tex = os.path.join(dir_name, f_base+".tex")
            file_pdf = os.path.join(dir_name, f_base+".pdf")
            tf = TemplateFiller(template="fc_tikz.tex")
            if self.fc_param["annotation.size"] is not None:
                fs = self.fc_param["annotation.size"]
            else:
                fs = 10         # Fallback to default
            figs = []           # List for variables pass to the jinja template
            for item in self.figure_list:
                # TODO: use relative path
                fig = item.subfig
                f_name = item.file_name
                sub_x_ = fig.get_x_fig()
                sub_y_ = fig.get_y_fig()
                sub_w_ = fig.get_w_fig()  # TODO: determine external image size?
                sub_h_ = fig.get_h_fig()
                bc = fig.fig_param["annotation.bc"]
                fc = fig.fig_param["annotation.fc"]
                alpha_ = fig.fig_param["annotation.alpha"]
                text_x_ = fig.get_x_text()
                text_y_ = fig.get_y_text()
                b_ = Bunch(
                    fig_x=fig.x_,
                    fig_y=fig.y_,
                    fig_w=fig.w_,
                    fig_h=fig.h_,
                    sub_x=sub_x_,
                    sub_y=sub_y_,
                    sub_w=sub_w_,
                    sub_h=sub_h_,
                    file_name=f_name,
                    alpha=alpha_,
                    text=fig.anno_text,
                    text_bc=bc,
                    text_fc=fc,
                    text_x=text_x_,
                    text_y=text_y_,
                    outline=outline,
                )
                figs.append(b_)

            tf.render(filename=file_tex,
                      font_size=fs,
                      font_name_main=self.fc_param["annotation.font"],
                      page_width=self.w,
                      page_height=self.h,
                      fig_list=figs)
            # TODO: use tikz engine to generate!
            # Possible to be running the subprocess with different TeX
            # directory. Must define the absolute path!
            subprocess.call(("latexmk",
                             "-xelatex",
                             "-pdf",
                             "-output-directory=%s" % dir_name,
                             file_tex))

default_figure_param = {"figure.lpad": 0.1,
                        "figure.rpad": 0.1,
                        "figure.tpad": 0.1,
                        "figure.bpad": 0.1,
                        "annotation.font": "Arial",
                        "annotation.size": None,
                        "annotation.weight": None,
                        "annotation.fc": None,
                        "annotation.bc": None,
                        "annotation.alpha": None,
                        "annotation.location": None}

valid_file_list = ("pdf",
                   "png",
                   "eps",
                   "jpg",
                   "jpeg")


class SubFigure(object):
    def __init__(self, loc, parent_param=None, fig_file=None, **argv):
        """
        `loc` is a 4-unit with (x_, y_, w_, h_)
        The SubFigure is better to be invoked by the FigureCollection,
        since the loc parameter is absolute value.
        """
        if len(loc) is not 4:
            raise ValueError("The loc takes a 4-unit item")
        try:
            self.x_ = float(loc[0])
            self.y_ = float(loc[1])
            self.w_ = float(loc[2])
            self.h_ = float(loc[3])
        except ValueError:
            warnings.warn("The loc parameter is invalid, abort!")
            raise
        self.fig_param = {}
        for key in default_figure_param:
            if key not in argv:
                if key in parent_param:
                    self.fig_param[key] = parent_param[key]
                else:
                    self.fig_param[key] = default_figure_param[key]
            else:
                self.fig_param[key] = argv[key]  # TODO: add assert!!
        # self.fig = fig
        self.anno_text = None
        self.fig_type = None

        """
        If user gives the `fig_file`, a file is used, 
        otherwise a mpl figure.
        Settings can be overwritten by `add_mpl_figure`
        or `add_file_figure`
        """
        if fig_file is not None:
            self.add_file_figure(fig_file)
        else:
            self.add_mpl_figure()

    def get_x_fig(self):
        return self.x_ + self.fig_param["figure.lpad"]*self.w_

    def get_y_fig(self):
        return self.y_ + self.fig_param["figure.tpad"]*self.h_

    def get_w_fig(self):
        return self.w_*(1 - self.fig_param["figure.lpad"]
                          - self.fig_param["figure.rpad"])

    def get_h_fig(self):
        return self.h_*(1 - self.fig_param["figure.tpad"]
                          - self.fig_param["figure.bpad"])

    def get_x_text(self):
        return self.x_ + self.fig_param["annotation.location"][0]*self.w_

    def get_y_text(self):
        return self.y_ + self.fig_param["annotation.location"][1]*self.h_

    def add_mpl_figure(self):
        """
        The size of the mpl figure is set only when plotting.
        Allowing user to change the geometry afterwards.
        """
        self.fig = Figure()
        self.backend = backend_agg.new_figure_manager_given_figure(1, self.fig)
        self.fig_type = "mpl"
        return self.fig

    def add_file_figure(self, filename):
        """
        Add a file instance as the figure.
        Supported file types are pdf, png, eps, jpeg, jpg,
        TODO: add tif or tiff support
        """
        if filename.strip().split(".")[-1] not in valid_file_list:
            raise TypeError("File type may not be supported!")
        self.fig = filename
        self.fig_type = "file"

    def set_plot_func(self, plot_func, *args, **argvs):
        """
        Register a plot_func to the fig instance.
        Plot is done only after pagesize is defined
        If tight_layout is used in plot_func,
        figure will give error after using set_size_inches!
        """
        if self.fig_type is "mpl" and self.fig is not None:
            self._plot_func = plot_func
            self._plot_params = Bunch(args=args, argvs=argvs)

    def save_fig(self, filename):
        """
        The geometry for each figure only effective when save it!
        """
        # self.w_fig = self.w_*(1 - self.fig_param["figure.lpad"]
                                # - self.fig_param["figure.rpad"])
        # self.h_fig = self.h_*(1 - self.fig_param["figure.tpad"]
                                # - self.fig_param["figure.bpad"])
        if self.fig_type is "mpl":
            self.fig.set_size_inches(self.get_w_fig(), self.get_h_fig())
            self._plot_fig()
            self.fig.savefig(filename)
        else:
            try:
                shutil.copyfile(self.fig, filename)
            except (shutil.SameFileError, FileNotFoundError, PermissionError):
                warnings.warn("The save file path is invalid!")
                raise

    def _plot_fig(self):
        """
        Use the user-defined function func to plot the mpl fig
        """
        args = self._plot_params.args
        argvs = self._plot_params.argvs
        if self.fig_type is "mpl" and self.fig is not None:
            self._plot_func(self.fig, *args, **argvs)

    """
    Some boring copy-paste functions for post setting
    """
    def set_lpad(self, lpad_):
        self.fig_param["figure.lpad"] = lpad_

    def set_rpad(self, rpad_):
        self.fig_param["figure.rpad"] = rpad_

    def set_tpad(self, tpad_):
        self.fig_param["figure.tpad"] = tpad_

    def set_bpad(self, bpad_):
        self.fig_param["figure.bpad"] = bpad_

    def set_pads(self, lpad=None, rpad=None, tpad=None, bpad=None):
        if lpad is not None:
            self.set_lpad(lpad)
        if rpad is not None:
            self.set_rpad(rpad)
        if tpad is not None:
            self.set_tpad(tpad)
        if bpad is not None:
            self.set_bpad(bpad)

    def set_geom(self, width=None, height=None):
        if width is not None:
            self.w_ = width
        if height is not None:
            self.h_ = height

    def set_annotation_text(self, s_):
        self.anno_text = s_

    def set_annotation_color(self, fc=None, bc=None, alpha=None):
        if fc is not None:
            self.fig_param["annotation.fc"] = fc
        if bc is not None:
            self.fig_param["annotation.bc"] = bc
        if alpha is not None:
            self.fig_param["annotation.alpha"] = alpha

    def set_annotation_loc(self, loc_):
        self.fig_param["annotation.location"] = loc_
