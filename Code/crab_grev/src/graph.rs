struct Node<'a> {
    name: String,
    neighbours: Vec<&'a str>,
}


struct Graph<'a> {
    nodes: Vec<Node<'a>>,
}

pub fn run() {
    let node1 = Node {
        name: String::from("node1"),
        neighbours: vec!["node1", "node2"],
    };

    let node2 = Node {
        name: String::from("node2"),

        neighbours: vec!["node1", "node2"],
    };

    let node3 = Node {
        name: String::from("node3"),
        neighbours: vec!["node1", "node2"],
        
    };

    let g1 = Graph {
        nodes: vec![node1, node2, node3]
    };
    for node in g1.nodes{
        println!("Node {} has the neighbours: {:?}", node.name, node.neighbours);
    }

}