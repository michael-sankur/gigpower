from . load_group import LoadGroup
from . capacitor import Capacitor


class CapacitorGroup(LoadGroup):
    dss_module_name = 'Capacitors'
    ele_class = Capacitor
