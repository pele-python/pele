For this project, as well as pele, we use the K&R c++ coding standard (more or
less) as a default

http://en.wikipedia.org/wiki/Indent_style

we use 4 spaces for indentation and *no tabs*

In generall, this means it looks like this ::

    int main(int argc, char *argv[])
    {
        ...
        while (x == y) {
            something();
            somethingelse();
     
            if (some_error) {
                do_correct();
            } else {
                continue_as_usual();
            }
        }
     
        finalthing();
        ...
    }

For a class it would be something like this::

    class MyClass{
    protected:
        double v1;
        double v2;
    public:
        MyClass()
            : v1(0),
              v2(0)
        { 
            do_something();
        }
    };

Some big things I noticed that we need to fix are

if and while statements have the brace on the same line::

    if (blah) {
        do something
    }

functions have braces on separate lines::

    void func()
    {
        do something
    }

always add white space after commas and around most operators (this is a pet peeve of mine ;)) )::

    func(var1,var2,var3);   // no!
    func(var1, var2, var3)  // yes!
    a=b+4; //no!
    a = b + 4; //yes!

with initializer lists, put the colon on a separate line from the name.  And the braces also::

    Minimizer::Minimizer()
        : val(0),
          ptr(1)
    {
        do something
    }

Try to keep the lines not too much longer than 80 characters.

Try to generally use braces with if statements.
for loops should always have braces.  it's just too dangerous otherwise.::

    // very easy to introduce problems
    if (condition)
        do_something;

add a space after for, if, while, etc::

    for(i = 0; i < N; ++i){ // no
    for (i = 0; i < N; ++i) { // yes

Put whitespace between operators::

    // way too little whitespace
    std::vector<energy_t> energies()const{return property_listing<energy_t>(&Minimum::energy);}
    // better, but it's still quite long and hard to read
    std::vector<energy_t> energies() const { return property_listing<energy_t>(&Minimum::energy); }
    // best
    std::vector<energy_t> energies() const 
    { 
        return property_listing<energy_t>(&Minimum::energy); 
    }


for naming, we follow the python convention and use CamelCase for class names
and lower_case_under_score for function and variable names.

Ideally functions should have a name that contains a verb. It should be the
action that the function performs.::

    value()     //bad
    get_value() //good

If you can't come up with a function name that describes what the function
does, that might be an indication that your function does too many things.
Try breaking up your function into many functions that each perform one action.
Functions should be simple enough that you can know what it's going to do just
by reading the name.

When naming the c++ tests we follow the standard convention and use ::

    TEST(ClassNameOrTestGroup, ActionPerformed_ExpectedResult)

the above is all CamelCase except for the single underscore separating the
action and the expected result.

For documentation, 
we try to follow the c++ Doxygen format.  That way we can
automatically generate nice looking documentation.  
http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html
In particular, functions should be documented like so::

    /**
     * return the energy
     */
    double get_energy();
