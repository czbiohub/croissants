#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 09:00:30 2019

@author: reganlamoureux
"""

class Example():
    classAttribute = 'Welcome'
    counter = 0
    
    def __init__(self, instanceAttribute, instanceAttribute2):
        self.instanceAttribute = instanceAttribute
        self.instanceAttribute2 = instanceAttribute2
        
    def setCounter(self, newCounter):
        self.counter = newCounter
        return self.counter
    
    def method(self, newSum):
        self.counter += newSum
        return self.counter
    
    def __str__(self):
        return "X is an object in class Example with attributes: {}, {}".format(self.instanceAttribute, self.instanceAttribute2)
    
    def __repr__(self):
        return "Example(" + str(self.instanceAttribute) + ',' + str(self.instanceAttribute2) +')'
    
class SubClass(Example):
    newClassAttribute = "Attribute for subclass"
    def __init__(self, instanceAttribute, instanceAttribute2, subClassAttribute):
        Example.__init__(self, instanceAttribute, instanceAttribute2)
        self.subClassAttribute = subClassAttribute
    
    def resetCounter(self):
        self.counter = 0
        return self.counter
    
    def __repr__(self):
        return "Example(" + str(self.instanceAttribute) + ',' + str(self.instanceAttribute2) + ',' +str(self.subClassAttribute)
        

class Dog():
    num_legs = 4
    counter = 0
    def __init__(self, breed, name, height):
        self.breed = breed
        self.name = name
        self.height = height
    def fetch(self):
        print("I am fetching")
        self.counter += 1
        return "Balls returned: {}".format(self.counter)
    def speak(self):
        print('woof')
    def __str__(self):
        return "Dog with breed {}, name {}, and height {}".format(self.breed, self.name, self.height)
    def __repr__(self):
        return "Dog({},{},{})".format(self.breed, self.name, self.height)
    
    
    
class AmericanCities():
    country = 'United States of America'
    def __init__(self, name, state, population):
        self.cityAttractions = []
        self.name = name
        self.state = state
        self.population = population
    
    def getName(self):
        return self.name
    
    def getState(self):
        return self.state
    
    def getPopulation(self):
        return self.population
    
    def increasePopulation(self, increase):
        self.population = self.population + increase
    
    def addCityAttractions(self, attraction):
        self.cityAttractions.append(attraction)
        
    def __repr__(self):
        return 'AmericanCities({},{},{},{})'.format(self.name, self.state, self.area, self.population)
    
class CaliforniaCities(AmericanCities):
    classCities = []
    def __init__(self, name, state, population, county, area):
        AmericanCities.__init__(self, name, state, population)
        self.area = area
        self.county = county
    def addClassCities(self):
        self.classCities.append(self.name)
        return self.classCities
    def predictWeather():
        print("No one can predict the weather in California")
    def __repr__(self):
        return "CaliforniaCities({},{},{},{},{})".format(self.name, self.state, self.population, self.county, self.area)
    
        
        
    
    
        
        