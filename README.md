## Network Abstractions

`network.py`: This module defines high-level abstractions for vertices in a network and for different kinds networks.

`gen_network.py`: This script generates a GraphML file representing a network having specified topology and size, and with the specified properties. See script for documentation on the allowed topologies and properties.

```bash
> python gen_network.py <topology> <size> <params> <outfile>
```

`graphml2gml.py`: This script converts a network file in GraphML format to one in GML format.

```bash
> python graphml2gml.py <infile>
```

## Software Dependencies

* [Python](https://www.python.org/)
* [Matplotlib](http://matplotlib.org/)
* [NetworkX](https://networkx.github.io/)
* [NumPy](http://www.numpy.org/)

## Contact

If you have any questions about the software, please email Swami Iyer at 
swami.iyer@gmail.com.
