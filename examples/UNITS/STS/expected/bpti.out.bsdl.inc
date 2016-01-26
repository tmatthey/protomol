//Auto generated include BSDL file by ProtoMol.
using 3D;
pointLight (0.8, [1,1,1]) {position [10000,5000,10000];};
ambientLight (0.3, [1,1,1]);
define b whitted {diffuse [1.0, 1.0, 1.0];ambient [0.5, 0.5, 0.5];};
define d whitted {diffuse [0.7, 0.7, 0.7];ambient [0.5, 0.5, 0.5];};
define c0 whitted {diffuse [1.0, 1.0, 1.0];ambient [0.5, 0.5, 0.5];};
define c1 whitted {diffuse [1.0, 0.1, 0.1];ambient [0.5, 0.1, 0.1];};
define c2 whitted {diffuse [0.1, 1.0, 0.1];ambient [0.1, 0.5, 0.1];};
define c3 whitted {diffuse [0.1, 0.1, 1.0];ambient [0.1, 0.1, 0.5];};
define c4 whitted {diffuse [0.1, 1.0, 1.0];ambient [0.1, 0.5, 0.5];};
define c5 whitted {diffuse [1.0, 1.0, 0.1];ambient [0.5, 0.5, 0.1];};
define c6 whitted {diffuse [1.0, 0.1, 1.0];ambient [0.5, 0.1, 0.5];};
define c7 whitted {diffuse [0.1, 0.1, 0.1];ambient [0.1, 0.1, 0.1];};
#define s sphere
#define c cylinder
#define u 0.25
#define v 0.1
grid [12,12,12]{
