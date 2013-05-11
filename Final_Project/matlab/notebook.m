M=csvread('/home/dk0r/git/comp-phys/Final_Project/midi.csv');
N=size(M,1);
P=2.2; %start delay

M(:,4) = 60 + round(40*rand(N,1));  % random volumes
M(:,5) = M(:,5);
%M(1,5) = 0;
M(:,6) = M(:,5) + round(40*rand(N,1)) - 1;
M(1,6) = 0; 
midi_new = matrix2midi(M);
writemidi(midi_new, 'testout2.mid');