
class: center, middle
# Python Classes for Beginners
## by: Regan Lamoureux
---

## Review: data types
* String: Character data; denoted with “”. (ex. `‘Hello world’` )
* Int: Integers (ex. 1)
* Float: Floating point numbers (ex 1.3)
* List: ordered sequence of objects (ex. `[1, 2, ‘hello’, world’]`)
* Boolean: Either True or False (ex. `‘World’.isalpha()`)

#### Some data types have built in functions
* `isalpha()` - checks if a string contains just letters (`‘hello’.isalpha()`)
* `pow(x, y)`-  returns x**y
* `.append(x)`- adds x to end of list
* `.remove(x)`- removes first instance of x in a list ([1,2,3].remove(3) *returns [1,2])

---

## Review: Functions
* Functions are a sequence of code that perform a specific task.
* Allows you to break up your code into more manageable pieces.
* Helpful with long code or if you are doing a specific task many times. 

```python
def function1(a): 
    X = 1+a
    return X

def function2(b):
    Y = 2+b
    return Y

...>>> value = function1(3) + function2(4)

...>>> value
       10
```

---


## What are Classes?
* Classes allows us to create new types or sets of objects.
* Attributes- variables or characteristics that you are able to obtain from your class
* Methods- operations on instances within the class
---

```python
class Example():
    classAttribute = ''
    counter = 0
    def __init__(self, instanceAttribute, instanceAttribute2):
        self.instanceAttribute = instanceAttribute
        self.instanceAttribute2 = instanceAttribute2
    def method(self, newSum):
        print('this is a class')
```

---

## Analogy

Think of it in the context of creating video game characters. Each character would have different attributes such as their name, age, where they are from, etc. Also each character could have different actions or powers. Some characters may be able to fly or cast spells. 

---

## Why are classes helpful?
* Keep related data together
* Keep track of objects with many attributes and methods
* Make code easier to read/look cleaner
---


```python
class Example():
    classAttribute = ''
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
        return "X is an object in class Example with attributes: {}, {}"
                .format(self.instanceAttribute, self.instanceAttribute2)
    
    def __repr__(self):
        return "Example(" + str(self.instanceAttribute) + ',' + 
                            str(self.instanceAttribute2) +')'
```
---

## Instance Attribute vs Class Attribute
* Instance Attributes are initialized with init or other methods and refer only to that specific object
* Class objects are initialized when an object is put into a class and is given to every object in that class.
* Be careful when assigning class and instance attributes

---

In the previous example:

   * classAttribute and counter are class attributes
    
   * instanceAttribute and instanceAttribute2 are instance attributes
    
If we initialized two objects in the example class:

    `x = Example('Hello', 'World')`
    
    `y = Example('Cupcakes', 'coding')`
    
What would be the class and instance attribute for each object?

---

## __init__
* `__init__` is the constructor method for python
* When you create an object within your class, `__init__` is automatically called.
* Most of the time used to set initial attributes.
* Self is a keyword that allows you to access a particular instance of the class.  


```python
def __init__(self, attribute1, attribute2):
    self.attribute1 = attribute1
    self.attriubte2 = attribute2
```
---

## __repr__
* `__repr__` is used to show string representation of a class
* 'Formal' notation
* The output should be what you would type in to get a result

```python
def __repr__(self):
    return "Example({},{})".format(self.instanceAttribute, 
                                   self.instanceAttribute2)
```
This should return `Example(x, y)` where x and y represent the two attributes initialized

---

## __str__
* very similar to `__repr__`
* `__str__`: ‘informal’ representation
* The output should be something easily readable

```python
def __str__(self):
        return "X is an object in class Example with attributes: {}, {}"
                .format(self.instanceAttribute, self.instanceAttribute2)
```
   
This should return `"X is an object in class Example with attributes: x, y"` where x and y are the two attributes initialized

---

# Let's Code a simple Class

We want to create a Dog class where the objects in this class represent dogs with:
* class attributes num_legs as 4 and counter as 0
* the characteristics of breed, name and height.
* a method called fetch where we print a statement and increase the counter to how many times this dog has fetched
* a speak method that prints out a bark
---

If you want, you can create a class called American Cities with:
* a class attribute of country as United States of America
* initialize the instance attributes of name, state, population and cityAttractions as an empty list
* create method getName which simply returns the name of that object
* create method getState which returns the state of that object
* create method getPopulation which returns the population of that object
* create method increasePopulation which takes an integer as an argument and adds that amount to the population of the object
* create a method that takes as an argument a string of an attraction in that city and appends that to cityAttractions
* Don't forget to add a `__repr__` and/or `__str__`
---

## Inheritance
* If you have multiple classes and certain properties overlap, you can create a class that is a modified version of an existing class.
* Syntax: `class SubClass(SuperClass):`
* The subclass inherits all the methods from the superclass
* With this you can add new methods or redefine methods from the superclass
---

```python
class SubClass(Example):
    newClassAttribute = "Attribute for subclass"
    def __init__(self, instanceAttribute, instanceAttribute2, 
                 subClassAttribute):
        Example.__init__(self, instanceAttribute, instanceAttribute2)
        self.subClassAttribute = subClassAttribute
    
    def resetCounter(self):
        self.counter = 0
        return self.counter
    
    def __repr__(self):
        return "Example({},{},{})".format(self.instanceAttribute, 
                                          str(self.instanceAttribute2, 
                                          str(self.subClassAttribute)
```
---

Lets Create a subclass of AmericanCities called CaliforniaCities with:
* a class attribute of an empty list 
* initialize AmericanCities with the instance attributes name, state, population, county, and area
* create a method called addClassCities. This will append the name of that object to the class list. This will create a complete list of all objects in that class.
* create a method called predictWeather that prints out the statement "No one can predict the weather in California"
* Don't forget your `__repr__`

```python
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
        return "CaliforniaCities({},{},{},{},{})".format(self.name, 
                self.state, self.population, self.county, self.area)
```


```python

```


```python

```
