# Fingerprint comparison algorithm based on graph theory

This fingerprint comparison algorithm defined on the concept of local structure is based on the article of Ratha N. K., Pandit V.D., Bolle R. M., Vaish V., “Robust Fingerprint Authentication Using Local Structural Similarity,” in *Proc. Workshop on Applications of Computer Vision*, pg. 29-34, 2000.
The fingerprint is represented by an undirected graph of the adjacent minutiae (Minutiae Adjacency Graph, MAG).
Given a fingerprint image in greyscale, it is assumed that the distinctive characteristics, in particular terminations and bifurcations, have already been extracted. The extracted minutiae are associated with the nodes of the graph, and the connections between the minutiae, if included in the neighbourhood considered to represent the arches.
The algorithm is divided into three phases, each of them provides a cost that is a parameter that indicates how much the nodes are coupled. This cost, combined with the other, allows determining the final score indicating the degree of similarity between the two fingerprints compared.
In the first phase, it is analysed the minutiae of a single fingerprint by selecting those which have similar characteristics in a neighbourhood. These will then be associated to those of the fingerprint, obtaining a set of pairs of nodes having similar characteristics.
In the second phase, it is checked whether the coupled nodes have similar characteristics, by comparing the distance between the respective nodes of each fingerprint.
In the third phase, it is extended the neighbourhood including the nodes excluded from the first phase, which they form arches with nodes of the first phase, which will be compared between the two fingerprints.
The algorithm considers the rotation, the translation, the partial overlapping and distortion of an image.

## Usage

The main.m file is working with [MATLAB](https://www.mathworks.com/products/matlab.html) (7.5 or more) and [Octave](https://www.gnu.org/software/octave/) (4.2 or more) with the [statistics](https://octave.sourceforge.io/) library.

## Author

* **Pietro Russo**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
