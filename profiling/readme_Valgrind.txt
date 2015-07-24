
Valgrind is a set of tool for examining what happens when code runs and profiling where bottlenecks lie. I'm familiar with two features of it: callgrind and memcheck, but there are others. It does this by simulating a CPU for the code to execute on. This makes it very, very, very slow. It also needs to be run on a compiled program: no ROOT macros in the interpreter. Fortunately, it's fairly easy to compile many ROOT macros*.

Code doesn't need to be compiled with any special options, although debug info (-g) can be a good idea. Turning off optimisation (-O0) can make things a little clearer, but much less useful. I always turn have optimisations on.

callgrind is the tool you use to see where the time is being spent in your code. We invoke it like this:
valgrind --tool=callgrind -v myprogram myargs

This then executes your code, apparently as normal but incredibly slowly, and saves information into a file named callgrind.out.pid where pid is the process ID. This process can be interupted at any time with ctrl-C, and you get all the information about the program's run up to that point. If you interupt with ctrl-Z, you'll lose the information, so don't do that.

The output file of callgrind is not really intended to be human readable - you need a decent GUI. I use an excellent Linux one called 'kcachegrind'. I understand that 'qcachegrind' is a similar program for Macs, but I've never used it.

kcachegrind is mostly fairly clear - you get a list where your code is spending its time, and what the breakdown of work in each function is. Work in a function is divided into other functions it calls, and internal work, called 'self'. Work that is considered 'self' is very basic things like operators on basic types, conditions, loops and memory management. The 'callee map' gives graphical breakdown of this. Each function also tells you how many times it was called, and from where. Occasionally, function names don't resolve (you just get hex addresses) but it still tells you what library its in and what called it, and you can usually work out what they are from context. I've only had this with the most basic of C standard functions (like exp).


memcheck is a tool for looking at memory use in programs. I mostly use it for finding memory leaks when necessary. I invoke with this:
valgrind --tool=memcheck --leak-check=yes -v myprog myargs > memcheck.txt
The output of this is human readable, so long as you use your text editor's find feature. 'definitely' is the most useful string to seach for, as in memory 'definitely lost'. Code compiled with debug flag is pretty much essential here, since it tells you where the leaked memory is assigned (and thus makes it pretty easy to find where it should have been deleted). Probably unneccessary for VaLOR for now.




* To compile ROOT macros to executables, you basically need three things:
1) A ROOT macro that obeys C++ rules, which they mostly should anyway.
2) int main()
3) the include statements.

The essential thing here is:
#ifndef __CINT__
This lets you have code that the compiler sees but regular ROOT execution ignores. If you encase all your include statements and your int main in these directives, you can use the same file as both a live macro and compilable source. If you're just profiling, you can usually call a macro with default arguments, so it's usually enough to add to the end of your file:
#ifndef __CINT__
int main(){
  mymacro();
  return 0;
}
#endif

