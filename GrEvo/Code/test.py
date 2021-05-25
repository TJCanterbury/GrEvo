import gressure as gr
import cProfile
print(gr.greet())
def reflect_n(node):
    """ returns the reflected node label """

    if node[0] == 'L':
        # If on left side add copy on right side
        node = 'R' + node[1:]
    
    elif node[0] == 'R':
        # If on right side add copy on left side
        node = 'L' + node[1:] 
    
    print(node)
numbers = list(range(0, 100000))
list = ["L" + str(e) for e in numbers]
cProfile.run("gr.print_vector(list)")
cProfile.run("for i in list: reflect_n(i)")
