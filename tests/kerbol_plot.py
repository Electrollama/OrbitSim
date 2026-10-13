__test__ = False

from UI.orbit_plotter import plot_system, plot_show
from Load.load import load_system

if __name__ == "__main__":
    plot_system(load_system('../Load/System Tables/Bodies_KSP.csv'), time=0.0)
    plot_show()