function w = li98cfs_wpot(t, d)
% potential function from li:98:cfs
t = t ./ d;
w = atan(t);
w(t == 0) = 1;
t(t == 0) = 1;
w = w ./ t;
