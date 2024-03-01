import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
from option_types import Option
from binomial_tree_options_pricing import binomial_tree

qtCreatorFile = "design.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self):        
        super().__init__()
        self.setupUi(self)
        
        self.calculate.clicked.connect(self.calculate)
        
    def calculate(self):

        if self.opt_type.text().upper() == "CALL": 
            opt_type = Option.CALL
        elif self.opt_type.text().upper() == "PUT": 
            opt_type = Option.PUT

        S0 = float(self.S0.text())
        rf = float(self.rf.text())
        strike = float(self.strike.text())
        time = float(self.time.text())
        n = int(self.time_steps.text())
        up_factor = 1 + float(self.up_factor.text())

        self.output.setText(str(round(binomial_tree(K = strike, 
                                                    T = time, 
                                                    S0 = S0, 
                                                    r = rf, 
                                                    N = n, 
                                                    u = up_factor, 
                                                    d = 1/up_factor, 
                                                    opttype = opt_type), 2)))
        

if __name__ == '__main__':

    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())