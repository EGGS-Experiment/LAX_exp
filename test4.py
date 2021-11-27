from artiq.dashboard import moninj
from sipyco.asyncio_tools import atexit_register_coroutine
from qasync import QEventLoop
from PyQt5 import QtCore, QtGui, QtWidgets

app = QtWidgets.QApplication(["ARTIQ Dashboard"])
loop = QEventLoop(app)

d_ttl_dds = moninj.MonInj()
loop.run_until_complete(d_ttl_dds.start(':1', 3250))
atexit_register_coroutine(d_ttl_dds.stop)

d_ttl_dds.ttl_dock.show()
d_ttl_dds.dds_dock.show()

app.exec_()
