The format of the instances located at the folders `campinas-*`, are as follows:

The number of nodes

```sh
VERTICES : 81
```

The number of arcs

```sh
ARCOS_NOREQ : 196
```

The number of blocks

```sh
BLOQUES : 5
```

The list of arcs

```sh
LISTA_ARISTAS_REQ : 
(node_i_id, node_j_id) coste i_j_cost
...
```

Where 
- `node_i_id`: The origin node id;
- `node_j_id`: The destiny node id;
- `i_j_cost`: The arc distance in meters.

And the blocks list

```sh 
LISTA_BLOQUE : 
node_id1, node_id2, ..., node_idk, block_profit
...
``` 

Where 
- `node_idk`: The block `k`-th node;
- `block_profit`: The block profit.
