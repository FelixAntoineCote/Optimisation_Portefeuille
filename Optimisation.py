import pandas as pd
import numpy as np
import datetime
from statsmodels.stats.weightstats import DescrStatsW
from scipy.optimize import minimize, Bounds, LinearConstraint

import warnings

# Disable pandas warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pandas")


class Optimisation():
    
    def __init__(self, file, riskFreeRate):
        self.df = pd.read_excel(file, index_col=0) + 1
        self.riskFreeRate = riskFreeRate
        
        columns = ["Date de debut", "Date de fin", "Rendement Attendu", "Variance Attendue", "Ratio de Sharp", 
                   "Commodities Value", "Commodities Momentum", "Commodities Multi-style", "Commodities Carry",
                   "US Stock Selection Value", "US Stock Selection Momentum", "US Stock Selection Defensive", "US Stock Selection Multi-style",
                   "Intl Stock Selection Value", "Intl Stock Selection Momentum", "Intl Stock Selection Defensive", "Intl Stock Selection Multi-style",
                   "Fixed income Value", "Fixed income Momentum", "Fixed income Defensive", "Fixed income Multi-style", "Fixed income Carry"
                   ]
        self.results = pd.DataFrame(columns=columns)
        
    def NegativeSharpRatio(self, x, returns):
        meanReturns = returns.mean().values
        covarientMatrix = returns.cov().values
        
        stdDeviation = np.sqrt(np.matmul(np.matmul(x, covarientMatrix), x.T))
        excessReturns = np.matmul(meanReturns,x.T)-self.riskFreeRate
        
        
        return -excessReturns/stdDeviation
        
        
    def OptimizeSharpRatio(self, returns):
        
        n = returns.shape[1]
        
        # Some constraints to block short selling and have sum of weights equal to 1
        bounds = Bounds(0,1)
        constraints = LinearConstraint(np.ones(n), 1,1)
        
        # Initial weights
        initialWeights = np.ones(n)/n
        
        # Optimization (minimize the negative sharp ratio)
        buff = minimize(self.NegativeSharpRatio, x0=initialWeights, args=(returns,), bounds=bounds, constraints=constraints, tol = 10**-4)
        
        if not buff.success:
            print(f"Optimization failed: {buff.message}")
        
        # What we're looking for
        startDate = returns.index[0]
        endDate = returns.index[-1]
        meanReturns = returns.mean().values
        covarientMatrix = returns.cov().values
        optimizedWeight = buff.x
        expectedReturns = np.dot(meanReturns, optimizedWeight.T)
        expectedVariance = np.matmul(np.matmul(optimizedWeight, covarientMatrix), optimizedWeight.T)
        sharpRatio = -self.NegativeSharpRatio(optimizedWeight, returns)
        
        return startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio
        
    def SimulateSharpRatio(self, window, markets, frequency = 3):
        # Isolate markets
        returns = self.df[markets]

        # Metrics
        startDateList = []
        endDateList = []
        optimizedWeightList = []
        expectedReturnsList = []
        expectedVarianceList = []
        sharpRatioList = []
        
        # Optimise the portfolio at selected frequency, default is every 3 months
        for index in range(0, len(returns), frequency):
                
            if index == 0:
                continue
            
            # To account for indexes smaller than the window
            windowIndex = max(index - window, 0)
            
            periodReturns = returns.iloc[windowIndex:index]
            startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio = self.OptimizeSharpRatio(periodReturns)
        
            startDateList.append(startDate)
            endDateList.append(endDate)
            optimizedWeightList.append(optimizedWeight)
            expectedReturnsList.append(expectedReturns)
            expectedVarianceList.append(expectedVariance)
            sharpRatioList.append(sharpRatio)
            
        return startDateList, endDateList, optimizedWeightList, expectedReturnsList, expectedVarianceList, sharpRatioList
            
        
            
    
    def AddToResults(self, startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio, markets, portfolioProportion):
        
        indexStartDate = self.results[self.results["Date de debut"] == startDate].index
        indexEndDate = self.results[self.results["Date de fin"] == endDate].index
        if startDate in self.results["Date de debut"].values and endDate in self.results["Date de fin"].values and set(indexStartDate) & set(indexEndDate):
            index = list(set(indexStartDate) & set(indexEndDate))
            self.results["Rendement Attendu"][index] = self.results["Rendement Attendu"][index] + expectedReturns*portfolioProportion
            self.results["Variance Attendue"][index] = self.results["Variance Attendue"][index] + expectedVariance*portfolioProportion
            self.results["Ratio de Sharp"][index] = self.results["Ratio de Sharp"][index] + sharpRatio*portfolioProportion
            
            for indexMarket, market in enumerate(markets):
                self.results[market][index] = optimizedWeight[indexMarket]*portfolioProportion
            
        else:
            newrow = {
                "Date de debut" : [startDate],
                "Date de fin": [endDate],
                "Rendement Attendu": [expectedReturns*portfolioProportion],
                "Variance Attendue": [expectedVariance*portfolioProportion],
                "Ratio de Sharp" : [sharpRatio*portfolioProportion]
            }
            for indexMarket, market in enumerate(markets):
                newrow[market] = optimizedWeight[indexMarket]*portfolioProportion
                
            buff = pd.DataFrame(newrow)
            self.results = pd.concat([self.results, buff], ignore_index=True)
            
    def ExportResults(self, fileName):
        self.results.to_excel(fileName)
                
            

        
            
# Monthly:
riskFreeRate = 2/(100*12) + 1
x = Optimisation("Data.xlsx", riskFreeRate)

commodities = ["Commodities Value", "Commodities Momentum", "Commodities Multi-style", "Commodities Carry"]
startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio = x.SimulateSharpRatio(36, commodities)
portfolioProportion = 0.57
for i in range(len(startDate)):
    x.AddToResults(startDate[i], endDate[i], optimizedWeight[i], expectedReturns[i], expectedVariance[i], sharpRatio[i], commodities, portfolioProportion)
    
equity = ["US Stock Selection Value", "US Stock Selection Momentum", "US Stock Selection Defensive", "US Stock Selection Multi-style",
          "Intl Stock Selection Value", "Intl Stock Selection Momentum", "Intl Stock Selection Defensive", "Intl Stock Selection Multi-style"]
startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio = x.SimulateSharpRatio(36, equity)
portfolioProportion = 0.28
for i in range(len(startDate)):
    x.AddToResults(startDate[i], endDate[i], optimizedWeight[i], expectedReturns[i], expectedVariance[i], sharpRatio[i], equity, portfolioProportion)
    
fixedIncome = ["Fixed income Value", "Fixed income Momentum", "Fixed income Defensive", "Fixed income Multi-style", "Fixed income Carry"]
startDate, endDate, optimizedWeight, expectedReturns, expectedVariance, sharpRatio = x.SimulateSharpRatio(36, fixedIncome)
portfolioProportion = 0.15
for i in range(len(startDate)):
    x.AddToResults(startDate[i], endDate[i], optimizedWeight[i], expectedReturns[i], expectedVariance[i], sharpRatio[i], fixedIncome, portfolioProportion)

x.ExportResults("Résultats_Partie1.xlsx")
