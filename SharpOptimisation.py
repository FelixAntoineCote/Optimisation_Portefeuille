import pandas as pd
from Optimisation import Optimisation




# Contraintes:
# Equity: 0.280 , Fixed Income: 0.1500 et Commodities: 0.5700
# Fenetre: 38 mois

# Import Data
file = "C:/Programation/Optimisation_Portefeuille/Data.xlsx"
df = pd.read_excel(file, index_col=0)

# Seperate the strategies
equity = ["US Stock Selection Value", "US Stock Selection Momentum", "US Stock Selection Defensive", "US Stock Selection Multi-style",
          "Intl Stock Selection Value", "Intl Stock Selection Momentum", "Intl Stock Selection Defensive", "Intl Stock Selection Multi-style"]
fixedIncome = ["Fixed income Value", "Fixed income Momentum", "Fixed income Defensive", "Fixed income Multi-style", "Fixed income Carry"]
commodities = ["Commodities Value", "Commodities Momentum", "Commodities Multi-style", "Commodities Carry"]

print(df[equity])