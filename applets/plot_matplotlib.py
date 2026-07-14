import numpy as np

from PyQt5.QtCore import Qt
import matplotlib
import matplotlib.pyplot as plt

from sipyco import pyon
from artiq.applets.simple import TitleApplet
from LAX_exp.applets.widget import QMainWindow

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt5.QtWidgets import QVBoxLayout, QWidget, QApplication
from sipyco.pc_rpc import Client
from PyQt5.QtCore import Qt, QTimer

class MatplotlibPlot(QMainWindow):
    """
    Class for creating an interactive Matplotlib plot in an applet

    Attributes:
        sc (MatplotlibCanvas): a matplotlib canvas to hold the plots
        data_points_list (list): list of data points clicked
        boxes (list): list of text boxes on figure
        box_clicked (ax.text): latest text box clicked
        default_keys (list): keys that will be used in the results dictionary passed to this class
    """

    def __init__(self, args, req, **kwargs):
        # Call super
        super().__init__(args, req, **kwargs)
        # create matplotlib canvas to insert plots into
        self.figure = Figure(figsize=(6, 8), constrained_layout=True)
        self.sc = FigureCanvasQTAgg(self.figure)

        self.setCentralWidget(self.sc)

        # instantiate clients to talk to dataset manager and master dashboard
        self.dataset_key = args.results
        self.dataset_db = Client(
            args.server,
            args.port_control,
            target_name = "dataset_db")
        self.ccb = Client(
            args.server,
            args.port_control,
            target_name="master_management"
        )

        self.applet_name = args.applet_name
        self.applet_group = args.applet_group

        # instantiate lists of boxes and points for later events
        self.data_points_list = []
        self.fit_lines = []
        self.boxes = []
        self.box_clicked = None

        # connect canvas to button press events
        self.sc.figure.canvas.mpl_connect('button_press_event', self.mouse_event)
        if args.title is not None:
            self.sc.figure.suptitle(args.title, fontsize=18)

        # default keys for dictionary passed to applet
        self.default_keys = ['x', 'y', 'z', 'errors', 'fit_x', 'fit_y', 'fit_z', 'subplot_x_labels', 'subplot_y_labels',
                             'subplot_titles', 'ylims', 'rid', 'legend_labels']

        if args.projection_3d is None:
            projection_3d = False
        elif args.projection_3d == "True":
            projection_3d = True
        else:
            projection_3d = False

        self.axes = self.create_axes(args.num_subplots, projection_3d)

    def create_axes(self, nsubplots=1, projection_3d = False):
        rows = nsubplots // 2
        rows = max(rows, 1)
        cols = nsubplots // rows + nsubplots % rows
        cols = max(cols, 1)

        if not projection_3d:
            axes = self.sc.figure.subplots(rows, cols)
        else:
            axes = self.sc.figure.subplots(rows, cols, subplot_kw = {'projection': "3d"})
        if not isinstance(axes, list) and not isinstance(axes, np.ndarray):
            axes = np.array([axes])
        return axes

    def update_applet(self, args):

        """
        Process run everytime dataset is modified and produces and updates the plot

        Args:
            args (dict): Dictionary of arguments
        """
        # clear old data
        for ax in self.axes.flatten():
            ax.cla()

        # # connect canvas to button press events
        # self.sc.figure.canvas.mpl_connect('button_press_event', self.mouse_event)

        # extract data from dictionary
        results = pyon.decode(self.get_dataset(args.results))

        # parse dictionary
        ys = self.get_from_dict(results, 'y', None)
        xs = self.get_from_dict(results, 'x', None)
        zs = self.get_from_dict(results, 'z', None)
        errors = self.get_from_dict(results, 'errors', None)
        fit_xs = self.get_from_dict(results, 'fit_x', None)
        fit_ys = self.get_from_dict(results, 'fit_y', None)
        fit_zs = self.get_from_dict(results, 'fit_z', None)
        x_labels = self.get_from_dict(results, 'subplot_x_labels', None)
        y_labels = self.get_from_dict(results, 'subplot_y_labels', None)
        z_labels = self.get_from_dict(results, 'subplot_z_labels', None)
        titles = self.get_from_dict(results, 'subplot_titles', None)
        ylims = self.get_from_dict(results, 'ylims', None)
        rid = self.get_from_dict(results, 'rid', None)
        legend_labels = self.get_from_dict(results, 'legend_labels', None)
        textbox_strs = self.get_from_dict(results,'textbox_strs', None)

        # determine number of datasets there are
        if zs is None:
            if len(ys.shape) == 1:
                num_datasets = 1
            else:
                num_datasets = ys.shape[0]
        else:
            if len(zs.shape) == 2:
                num_datasets = 1
            else:
                num_datasets = zs.shape[0]


        # ensure number of datasets is equal to number of subplots
        if args.num_subplots not in (1, num_datasets):
            raise ValueError(
                "num_subplots must be either 1 or equal to the number of datasets"
            )

        if args.num_subplots == 1:
            ax_ind = 0
            if zs is not None:
                self.plot(
                    xs,
                    ys,
                    errors,
                    fit_xs,
                    fit_ys,
                    z = zs,
                    fit_z = fit_zs,
                    x_label = x_labels,
                    y_label=y_labels,
                    title=titles,
                    ylim=ylims,
                    ind=ax_ind,
                    rid=rid,
                    legend_label=legend_labels,
                    textbox_str = textbox_strs[0]
                )

            elif len(ys.shape) == 1:
                self.plot(
                    xs,
                    ys,
                    errors,
                    fit_xs,
                    fit_ys,
                    x_label=x_labels,
                    y_label=y_labels,
                    title=titles,
                    ylim=ylims,
                    ind=ax_ind,
                    rid=rid,
                    legend_label=legend_labels,
                    textbox_str = textbox_strs[0]
                )
            else:
                for ind, y in enumerate(ys):
                    x = self.get_plot_element(xs, ind)
                    error = self.get_plot_element(errors, ind)
                    fit_x = self.get_plot_element(fit_xs, ind)
                    fit_y = self.get_plot_element(fit_ys, ind)
                    legend_label = self.get_plot_element(legend_labels, ind)
                    textbox_str = self.get_plot_element(textbox_strs, ind)

                    x_label = x_labels.item() if x_labels is not None and x_labels.shape == () else x_labels
                    y_label = y_labels.item() if y_labels is not None and y_labels.shape == () else y_labels
                    title = titles.item() if titles is not None and titles.shape == () else titles
                    ylim = ylims

                    self.plot(
                        x,
                        y,
                        error,
                        fit_x,
                        fit_y,
                        x_label = x_label,
                        y_label=y_label,
                        title=title,
                        ylim=ylim,
                        ind=ax_ind,
                        rid=rid,
                        legend_label=legend_label,
                        textbox_str = textbox_str
                    )

        else:
            for ind, y in enumerate(ys):
                x = self.get_plot_element(xs, ind)
                error = self.get_plot_element(errors, ind)
                fit_x = self.get_plot_element(fit_xs, ind)
                fit_y = self.get_plot_element(fit_ys, ind)
                x_label = self.get_plot_element(x_labels, ind)
                y_label = self.get_plot_element(y_labels, ind)
                title = self.get_plot_element(titles, ind)
                ylim = self.get_plot_element(ylims, ind)
                legend_label = self.get_plot_element(legend_labels, ind)
                textbox_str = self.get_plot_element(textbox_strs, ind)

                self.plot(
                    x,
                    y,
                    error,
                    fit_x,
                    fit_y,
                    x_label=x_label,
                    y_label =y_label,
                    title = title,
                    ylim=ylim,
                    ind=ind,
                    rid=rid,
                    legend_label=legend_label,
                    textbox_str = textbox_str
                )

        # set the matplotlib plot in the Qt widget and show the widget
        self.sc.figure.set_constrained_layout(True)
        self.sc.draw_idle()

    """PLOTTING FUNCTIONS"""

    def plot(self, x, y, error, fit_x, fit_y, x_label="", y_label="", z_label = "",
             title="", ylim = None, ind=0, rid=None,
             z=None, fit_z = None, textbox_str = None,
             legend_label = None):
        """
        Plot the data in a matplotlib subplots

        Args:
            x : x-data
            y : y-data
            error : error bars
            fit_x : domain for the fitted data
            fit_y : the fit for the data
            x_label : label for the x-axis of the subplot
            y_label : label for the y-axis of the subplot
            z_label: label for the z-axis of the subplot
            title : title for the subplot
            ind: index (location) of the subplot in the matplotlib figure
            rid: run number of experiment
            ylim: y limits
            textbox_str: text to write to figure
        """

        # verify the variables are of the correct type
        x = self.verify_aux_arg_type(y, x, "X")
        error = self.verify_aux_arg_type(y, error, "Error")
        fit_y = self.verify_aux_arg_type(fit_y, fit_y, "Fit Y")

        if fit_y is not None and fit_x is not None:
            if len(fit_y) > len(fit_x):
                fit_y = None
            elif len(fit_x) > len(fit_y):
                # trim fit x data
                fit_x = fit_x[:len(fit_y)]
        elif fit_y is not None and x is not None:
            fit_x = x
        else:
            fit_x = None

        # get features of the plot legend
        handles, labels = self.axes[ind].get_legend_handles_labels()

        # only plot if a new experiment is run
        if z is None:
            x_list = [x]
            y_list = [y]
            error_list = [error]
            for data_ind, xx in enumerate(x_list):
                data_points = self.axes[ind].errorbar(xx, y_list[data_ind], error_list[data_ind], marker="o", linestyle="-", label=legend_label,
                                                  markersize = 6)
                self.data_points_list.append(data_points)
        elif z is not None:
            x_list = [x]
            y_list = [y]
            z_list = [z]
            for data_ind, xx in enumerate(x_list):
                data_points = self.axes[ind].scatter(xx, y_list[data_ind], z_list[data_ind], marker="o", label=rid)
                self.data_points_list.append(data_points)


        # plot fit to data and set titles and labels
        if fit_y is not None and fit_x is not None and fit_z is None:
            fit_x_list = [fit_x]
            fit_y_list = [fit_y]
            for fit_ind, f_x in enumerate(fit_x_list):
                fit_lines = self.axes[ind].plot(f_x, fit_y_list[fit_ind])
                self.fit_lines.append(fit_lines)
        elif fit_y is not None and fit_x is not None and fit_z is not None:
            fit_x_list = [fit_x]
            fit_y_list = [fit_y]
            fit_z_list = [fit_z]
            for fit_ind, f_x in enumerate(fit_x_list):
                fit_lines = self.axes[ind].plot_wireframe(f_x, fit_y_list[fit_ind], fit_z_list[fit_ind])
                self.fit_lines.append(fit_lines)
        if title is not None:
            self.axes[ind].set_title(f"{title} \n "
                                     f"RID: {rid}")
        else:
            self.axes[ind].set_title(f"RID: {rid}")
        if x_label is not None and hasattr(self.axes[ind], "set_xlabel"):
            self.axes[ind].set_xlabel(x_label)
        if y_label is not None and hasattr(self.axes[ind], "set_ylabel"):
            self.axes[ind].set_ylabel(y_label)
        if z_label is not None and hasattr(self.axes[ind], "set_zlabel"):
            self.axes[ind].set_zlabel(z_label)
        if ylim is not None:
            self.axes[ind].set_ylim(np.min(ylim), np.max(ylim))
        if hasattr(self.axes[ind], 'text') and textbox_str is not None:
            # Define the text box properties
            box_properties = dict(
                boxstyle='round',
                facecolor='wheat',  # Inside color of the box
                edgecolor='black',  # Border color
                alpha=0.75,  # Transparency (0 = invisible, 1 = solid)
                pad=0.25  # Space between the text and the box border
            )
            self.axes[ind].text(0.05,0.95, textbox_str, fontsize = 16,transform=self.axes[ind].transAxes, va='top', ha='left',
                                bbox=box_properties)

        # style plot
        self.axes[ind].ticklabel_format(axis='x', style='plain', useOffset=False)
        self.axes[ind].ticklabel_format(axis='y', style='plain', useOffset=False)
        # set xticks by removing Nones from list
        x_filtered = [x_ele for x_ele in x if x_ele is not None]
        # self.axes[ind].set_xticks(np.round(np.linspace(np.min(x_filtered), np.max(x_filtered), 5),3))
        if legend_label is not None:
            self.axes[ind].legend()
        self.axes[ind].grid(True)

    """VERIFICATION AND HELPER FUNCTIONS"""

    def get_from_dict(self, d, key, default=None):
        if key in list(d.keys()):
            return np.array(d[key])
        else:
            return default

    def get_plot_element(self, plot_element, ind):
        """
        Get an element (title, xlabel, data, etc.) from a list or returns the element

        Args:
            plot_element : list of element (data or label)
            ind: index of the label in the list
        """
        # see if plot label comes in array
        if isinstance(plot_element, np.ndarray):
            # check if array is filled with all Nones
            if (plot_element == None).all():
                plot_element = None
            elif len(plot_element.shape) < 1:
                plot_element = plot_element
            # otherwise just grab the label
            elif ind > len(plot_element)-1:
                plot_element = None
            else:
                plot_element = plot_element[ind]
        elif isinstance(plot_element, list) and ind < len(plot_element)-1:
            plot_element = plot_element[ind]

        return plot_element

    def verify_aux_arg_type(self, y, aux_arg, aux_arg_name=""):
        """
        Verify the type of the not required datasets

        Args:
            y: the required dataset
            aux_arg: the not required dataset
            aux_arg_name: the name of the not required dataset
        """
        if aux_arg is not None:
            if (aux_arg == None).all():
                aux_arg = None

            elif not len(aux_arg):
                aux_arg = None
            else:
                if aux_arg.shape != y.shape:
                    self.logger.warning(f'{aux_arg_name} data does not match shape of Y data')
                    aux_arg = None

        return aux_arg

    """EVENT HANDLING"""

    def mouse_event(self, event):
        """
        Record a mouse click on data point or text box

        Args:
            event (Event): event recording a mouse click
        """

        # if click was on a text box already created then record and don't create new text box
        box_clicked = [box_clicked for box_clicked in self.boxes if box_clicked.contains(event)[0]]
        if len(box_clicked) > 0:
            self.box_clicked = box_clicked[0]
        # if data point was clicked create text box and record
        else:
            self.box_clicked = None
            artist_clicked = [picked.lines[0] for picked in self.data_points_list if picked.lines[0].contains(event)[0]]
            if len(artist_clicked) > 0:
                self.boxes.append(artist_clicked[0].axes.text(event.xdata, event.ydata,
                                                              f"x: {event.xdata:.4f} \ny: {event.ydata:.4f}",
                                                              fontsize=12,
                                                              bbox={'facecolor': 'white', 'pad': 4,
                                                                    'edgecolor': 'black',
                                                                    'alpha': 1}))

            # check fit lines
            else:
                self.box_clicked = None
                artist_clicked = [line for lines in self.fit_lines for line in lines if line.contains(event)[0]]

                if len(artist_clicked) > 0:
                    self.boxes.append(artist_clicked[0].axes.text(event.xdata, event.ydata,
                                                                  f"x: {event.xdata:.4f} \ny: {event.ydata:.4f}",
                                                                  fontsize=12,
                                                                  bbox={'facecolor': 'white', 'pad': 4,
                                                                        'edgecolor': 'black',
                                                                        'alpha': 1}))

        self.sc.draw_idle()

    def keyPressEvent(self, event):
        """
        If a text box was clicked and then the delete or backspace key is pressed, remove the text box

        Args:
            event (Event): An event object recording a key press
        """
        if event.key() == Qt.Key_Backspace or event.key() == Qt.Key_Delete:
            # remove clicked box from list of boxes after key event
            self.boxes = [box for box in self.boxes if box != self.box_clicked]
            # if box was clicked prior to key event remove it
            if self.box_clicked is not None:
                self.box_clicked.remove()
                self.box_clicked = None

        self.sc.draw_idle()


    def closeEvent(self, event):
        """
        Delete this applet's run-specific plotting dataset when the window closes.
        """
        dataset_key = getattr(self, "dataset_key", None)

        try:
            if dataset_key is not None and '.rid_' in dataset_key:
                # delete temporary dataset
                self.dataset_db.delete(dataset_key)
                self.dataset_db.close_rpc()
        except Exception as e:
            self.logger.warning("Dataset cleanup failed: %s", repr(e))


def main():
    # Create applet object
    applet = TitleApplet(MatplotlibPlot)

    # get arguments
    applet.argparser.add_argument("--x-label", type=str, default="", required=False)
    applet.argparser.add_argument("--y-label", type=str, default="", required=False)
    applet.argparser.add_argument("--num-subplots", type=int, default=1, required=False)
    applet.argparser.add_argument("--projection_3d", type=str, default=False, required=False)
    applet.argparser.add_argument("--applet-name", type=str, default=None)
    applet.argparser.add_argument(
        "--applet-group",
        nargs="*",
        default=None
    )

    applet.argparser.add_argument("--big-applet", action="store_true")
    applet.argparser.add_argument("--delete-on-close", action="store_true")

    # get datasets
    applet.add_dataset("results", "dictionary of experimental results")

    # run applet
    applet.run()

if __name__ == "__main__":
    main()
