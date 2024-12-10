import pandas as pd
import numpy as np
import datetime


class Optimisation():
    
    def __init__(self, file):
        self.df = pd.read_excel(file, index_col=0) + 1
        self.riskFreeRate = 1.02/12
        
    def SharpRatio(self, weights, date, window, markets):
        
        # Get the returns of interest
        indexDate = self.df.index.get_loc(date)
        returns = self.df[markets][indexDate-window:indexDate]

        averageReturns = returns.mean().values
        covarianceMatrix = returns.cov().values
        weights = np.transpose(weights)
        
        portfolioReturns = np.dot(averageReturns, weights)
        portfolioVariance = np.dot(weights.T, np.dot(covarianceMatrix, weights))
        portfolioStdDeviation = np.sqrt(portfolioVariance)
        
        SharpRatio = (portfolioReturns - self.riskFreeRate)/portfolioStdDeviation
        
        return SharpRatio
    
    def OptimiseSharpRatio(self, weights, date, window, markets, riskFreeRate = 1.02/12):
        
        
        
        
            
        
        
x = Optimisation("C:/Programation/Optimisation_Portefeuille/Data.xlsx")

commodities = ["Commodities Value", "Commodities Momentum", "Commodities Multi-style", "Commodities Carry"]
weights = np.ones(len(commodities))
date = "2022-10-31"

x.SharpRatio(weights, date, 30, commodities)
        
