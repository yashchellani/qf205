import sys
from PyQt5.QtWidgets import QMainWindow, QApplication
from PyQt5 import uic
import logging

from option_types import Option
from binomial_tree_options_pricing import binomial_tree
from monte_carlo_options_pricing import monte_carlo_option_price

from blackScholes import calculate_ftcs,calculate_crank_nicolson,calculate_exact_black_scholes

qtCreatorFile = "design.ui"
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)

# Setup logging
logger = logging.getLogger('OptionPricingApp')
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self):        
        super().__init__()
        self.setupUi(self)
        
        self.choose.clicked.connect(self.choose_method)

        self.calculate.clicked.connect(self.calculate_out)

    def choose_method(self):
        method = self.choice.currentText()
        self.S0.setText("S0")
        self.rf.setText("risk-free rate")
        self.strike.setText("strike price")
        self.time.setText("time to maturity (in years)")
        self.opt_type.setText("put/call")
        
        if method == "Binomial Tree":
            self.wildcard1.setText("time steps")
            self.wildcard2.setText("upmove factor (e.g. 0.1)")

        elif method == "Monte Carlo":
            self.wildcard1.setText("sigma aka volatility") 
            self.wildcard2.setText("no. of simulations")

        else:  # Applies to Crank-Nicolson, FTCS, and Black-Scholes
            self.wildcard1.setText("sigma aka volatility")
            self.wildcard2.hide()  # Hide second wildcard if not needed
        

    def calculate_out(self):
        try: 
            S0 = float(self.S0.text()) 
            rf = float(self.rf.text()) 
            strike = float(self.strike.text()) 
            time = float(self.time.text()) 
            method = self.choice.currentText() 
 
            if self.opt_type.text().upper() == "CALL":  
                opt_type = Option.CALL 
            elif self.opt_type.text().upper() == "PUT":  
                opt_type = Option.PUT 
 
            # Define a dictionary to act as a switch-case 
            method_function = { 
                "Binomial Tree": lambda: binomial_tree( 
                    K=strike, T=time, S0=S0, r=rf,  
                    N=int(self.wildcard1.text()),  
                    u=1 + float(self.wildcard2.text()),  
                    d=1 / (1 + float(self.wildcard2.text())),  
                    opttype=opt_type 
                ), 
                "Monte Carlo": lambda: monte_carlo_option_price( 
                    S0=S0, K=strike, T=time, r=rf,  
                    sigma=float(self.wildcard1.text()),  
                    num_simulations=int(self.wildcard2.text()),  
                    opttype=opt_type 
                ), 
                "Crank-Nicolson": lambda: calculate_crank_nicolson( 
                    S0, strike, time, rf, float(self.wildcard1.text()) 
                ), 
                "FTCS": lambda: calculate_ftcs( 
                    S0, strike, time, rf, float(self.wildcard1.text()) 
                ), 
                "Black-Scholes": lambda: calculate_exact_black_scholes( 
                    S0, strike, time, rf, float(self.wildcard1.text()) 
                ) 
            } 
 
            # Execute the appropriate function based on the user's choice 
            result = method_function[method]()  
            logger.info('result: ' + str(result)) 
            self.output.setText(str(round(result, 4)))

        except Exception as e:
            logger.error(f"Error calculating option price: {e}")
                

if __name__ == '__main__':

    app = QApplication(sys.argv)
    main = Main()
    main.show()
    sys.exit(app.exec_())