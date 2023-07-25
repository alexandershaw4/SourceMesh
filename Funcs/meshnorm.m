function N = meshnorm(g)

tr = triangulation(double(g.faces),double(g.vertices(:,1)),...
                                   double(g.vertices(:,2)),...
                                   double(g.vertices(:,3)) );
g.FaceNorm   = tr.faceNormal;
N            = -double(tr.vertexNormal);           

normN = sqrt(sum(N.^2,2));
normN(normN < eps) = 1;
N     = N ./ repmat(normN,1,3);
g.VertexNorm = N;

end