try:
    from viper.viper_dolfin import *
except:
    def plot(*args):
        raise RuntimeError, "Unable to plot (Viper plotter not available)"

    class Plotter:
        def __getattr__(self, name):

            def method(*args, **kwargs):
                raise RuntimeError, "Unable to plot (Viper plotter not available)"

            return method
