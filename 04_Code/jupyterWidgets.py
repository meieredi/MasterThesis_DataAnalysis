import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display, Javascript


def runAllWidget():
    def run_all(ev):
        display(Javascript('IPython.notebook.execute_cell_range(IPython.notebook.get_selected_index()+2, IPython.notebook.ncells())'))
    
    button = widgets.Button(description="Run all cells below")
    button.on_click(run_all)

    print('Finally, press "Run all cells below" button below!')
    
    return button
    

def plotWidget():
    widget_plotDLS = widgets.Checkbox(
                    description='Plot DLS', 
                    value=True)
    widget_plotMDI = widgets.Checkbox(
                    description='Plot MDI',
                    value=True)
    widget_onlyOverview = widgets.Checkbox(
                    description='Only overview plots',
                    value=True)

    print("Please decide what results shall be plotted:")

    box = widgets.VBox([widget_plotDLS, widget_plotMDI, widget_onlyOverview])

    return [box, widget_plotDLS, widget_plotMDI, widget_onlyOverview]


def defaultWidgets():

    # Define relevant widgets

    checkbox_showTitle = widgets.Checkbox(
                        description='Show title on plots', 
                        value=True)

    menu_cmap = widgets.Dropdown(
                options=['coolwarm', 'copper', 'viridis', 'plasma', 'inferno', 'cividis'],
                value='coolwarm',
                description='Colormap:')

    checkbox_cmap = widgets.Checkbox(
                description='Invert colormap', 
                value=True)

    menu_font = widgets.Dropdown(
                options=['Arial', 'LaTeX', 'Default'],
                value='Arial',
                description='Font:')

    slider_fontSize = widgets.FloatSlider(
                    value=15,
                    min=1,
                    max=20.0,
                    step=1,
                    description='Fontsize:')

    print("Please select the settings for plots:")
    print("(Recommended font size is 12 for documents and 18 for presentations)\n") 

    box = widgets.VBox([checkbox_showTitle, menu_cmap, checkbox_cmap, menu_font, slider_fontSize])

    return [box, checkbox_showTitle, menu_cmap, checkbox_cmap, menu_font, slider_fontSize]


def setDefault(onlyOverview, showTitle, cmap, invert, font, fontSize):
    
    # Old version (w/o widgets):
    # plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.coolwarm(np.linspace(0,1,4)))

    # New version (w/ widgets)
    if invert:
        cmap += '_r'

    if font == 'LaTeX':
        plt.rcParams.update({
        'text.usetex': True,
        'font.family': 'serif',
        'font.size'  : fontSize,
        'font.serif': ['Computern Modern Roman'],
        })
    elif font == 'Arial':
        plt.rcParams.update({
        'text.usetex': False,
        'font.family': 'sans-serif',
        'font.size'  : fontSize,
        'font.sans-serif': ['Arial'],
        })
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.rm'] = 'Arial'
        plt.rcParams['mathtext.it'] = 'Arial'
    else:
        plt.style.use('default')
        plt.rcParams.update({
        'font.size'  : fontSize
        })
    if fontSize >= 15:
        plt.rcParams.update({
        'figure.titlesize'  : 20
        })
    else:
        plt.rcParams.update({
        'figure.titlesize'  : 16
        })
        
    print('Plot:       overview') if onlyOverview else print('Plot:       everything')
    print('Title:     ', showTitle)
    print('Colormap:  ', cmap)
    print('Font:      ', font)
    print('Fontsize:   {}'.format(fontSize))

    return cmap

def pathWidgets():

    text_meiered = widgets.Text(
                   value='Data/meiered/DLS',
                   description='Path meiered')
    text_gerltm = widgets.Text(
                   value='Data/gerltm/DLS',
                   description='Path gerltm')
    
    box = widgets.VBox([text_meiered, text_gerltm])

    return [box, text_meiered.value, text_gerltm.value]