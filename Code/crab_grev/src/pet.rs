use petgraph::graph::{UnGraph};

pub fn run() {


    // Create an undirected graph with `i32` nodes and edges with `()` associated data.
    let g = UnGraph::<i32, ()>::from_edges(&[
        (1001, 1102), (1102, 1003), (1202, 1001), 
        (1202, 1003)]);
    let n1 = g.neighbors(1001.into());
    for n in n1{
        println!("{:?}", n);
    }
}