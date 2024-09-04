## Campinas

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

## Limoeiro & Alto Santo

The format of the instance located at the folder `carlos`, are as follows:

In the first line we have three numbers

```sh
#n #a #b
```

where
- `#n`: The number of nodes;
- `#a`: The number of arcs;
- `#b`: The number of blocks.

The next `#n` lines, are the graph nodes

```sh
N id1 x1 y1 -1
N id2 x2 y2 b1,...,bk
```

where
- `id`: Stands for the node id;
- `x`: Is the node x-axis value;
- `y`: Is the node y-axis value;
- `bi`: Is the `i`-th block id containing node `id`. Case no blocks containing the node, the entry `-1` is used.

The next `#a` lines, are the graph arcs

```sh
A orig_id dest_id dist block_id #cases xxxxx 
```

where
- `orig_id`: Is the arc origin node;
- `dest_id`: Is the arc destiny node;
- `dist`: Is the arc distance;
- `block_id`: Is the block id containing the arc, case the arc does not belong to any block, you will find -1;
- `#cases`: Is the number of reported cases;
