function d=opSmoothUnion(d1, d2, k )
    h = max( k-abs(d1-d2), 0.0 );
    d = min( d1, d2 ) - h.*h.*0.25./k;
end