"""
moninj: try taking docks
"""
from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QMainWindow, QGridLayout
from artiq.dashboard import moninj
from artiq.frontend.artiq_dashboard import MainWindow
from qasync import QEventLoop
import asyncio
import time

app = QApplication(["scde title"])
loop = QEventLoop(app)
asyncio.set_event_loop(loop)

th1=moninj.MonInj()
loop.run_until_complete(th1.start('::1', 3250))
ttl1=th1.ttl_dock

mw = MainWindow('::1')
mw.addDockWidget(QtCore.Qt.RightDockWidgetArea, ttl1)
mw.show()

loop.run_until_complete(mw.exit_request.wait())