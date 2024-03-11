import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic

from option_types import Option
from binomial_tree_options_pricing import binomial_tree
from monte_carlo_options_pricing import monte_carlo_option_price

qtCreatorFile = "design.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self):        
        super().__init__()
        self.setupUi(self)
        
        self.choose.clicked.connect(self.choose_method)

        self.calculate.clicked.connect(self.calculate_out)

    def choose_method(self):
        if self.choice.currentText() == "Binomial Tree":
            self.S0.setText("S0")
            self.rf.setText("risk-free rate")
            self.strike.setText("strike price")
            self.time.setText("time to maturity (in years)")
            self.opt_type.setText("put/call")

            # use wildcard variables because the methods dont take in the exact same parameters
            self.wildcard1.setText("time steps")
            self.wildcard2.setText("upmove factor (e.g. 0.1)")


        elif self.choice.currentText() == "Monte Carlo":
            self.S0.setText("S0")
            self.rf.setText("risk-free rate")
            self.strike.setText("strike price")
            self.time.setText("time to maturity (in years)")
            self.opt_type.setText("put/call")

            self.wildcard1.setText("sigma aka volatility") 
            self.wildcard2.setText("no. of simulations") 
            

        elif self.choice.currentText() == "Crank-Nicolson":
            pass
        

    def calculate_out(self):

        if self.opt_type.text().upper() == "CALL": 
            opt_type = Option.CALL
        elif self.opt_type.text().upper() == "PUT": 
            opt_type = Option.PUT


        if self.choice.currentText() == "Binomial Tree":
            S0 = float(self.S0.text())
            rf = float(self.rf.text())
            strike = float(self.strike.text())
            time = float(self.time.text())

            n = int(self.wildcard1.text())
            up_factor = 1 + float(self.wildcard2.text())

            self.output.setText(str(round(binomial_tree(K = strike, T = time, S0 = S0, 
                                                        r = rf, N = n, u = up_factor, 
                                                        d = 1/up_factor, opttype = opt_type), 2)))
            

        if self.choice.currentText() == "Monte Carlo":
            S0 = float(self.S0.text())
            rf = float(self.rf.text())
            strike = float(self.strike.text())
            time = float(self.time.text())

            sigma = float(self.wildcard1.text())
            num_simulations = 1 + int(self.wildcard2.text())

            self.output.setText(str(round(monte_carlo_option_price(S0 = S0, K = strike, T = time, r = rf, 
                                                                   sigma = sigma, num_simulations = num_simulations, 
                                                                   opt_type = opt_type), 2)))
        

if __name__ == '__main__':

    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())