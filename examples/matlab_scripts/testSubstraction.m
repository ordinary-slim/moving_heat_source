bigMesh = Mesh( -1, 1, 10 );
smallMesh = Mesh( 0.0, 1, 5 );

bigMesh.activeElements
bigMesh.substract( smallMesh );
bigMesh.activeElements