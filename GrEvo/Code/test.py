import gressure as gr
import cProfile
numbers = list(range(0, 10))
list = ["L" + str(e) for e in numbers]
print(gr.Reflect(list))