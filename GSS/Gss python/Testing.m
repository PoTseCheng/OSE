Qn =10; % The default option
eps = [3.436159118837738;2.532731674232790;1.756683649299882;1.036610829789514;0.3429013272237046;-0.3429013272237046;-1.036610829789514;-1.756683649299882;-2.532731674232790;-3.436159118837738];
weight = [7.640432855232621e-06;0.001343645746781233;0.03387439445548106;0.2401386110823147;0.6108626337353258;0.6108626337353258;0.2401386110823147;0.03387439445548106;0.001343645746781233;7.640432855232621e-06];

N=3

n_nodes = Qn^N;        % Total number of integration nodes (in N dimensions)

z1 = zeros(n_nodes,N); % A supplementary matrix for integration nodes; 
                       % n_nodes-by-N 
w1 = ones(n_nodes,1);  % A supplementary matrix for integration weights; 
                       % n_nodes-by-1

for i = 1:N            
   z1i = [];           % A column for variable i to be filled in with nodes 
   w1i = [];           % A column for variable i to be filled in with weights 
   for j = 1:Qn^(N-i)
       for u=1:Qn
           z1i = [z1i;ones(Qn^(i-1),1)*eps(u)];
           w1i = [w1i;ones(Qn^(i-1),1)*weight(u)];
       end
   end
   z1(:,i) = z1i;      % z1 has its i-th column equal to z1i 
   w1 = w1.*w1i;       % w1 is a product of weights w1i 
end

z = sqrt(2).*z1;       % Integration nodes; n_nodes-by-N; for example, 
                       % for N = 2 and Qn=2, z = [1 1; -1 1; 1 -1; -1 -1]

w = w1/sqrt(pi)^N;