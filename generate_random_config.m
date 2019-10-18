
rng(1);
num_samples_to_generate = 20;
min_val = 0;
max_val = 2*pi;
num_DOF = 5;
start_config = [1,num_DOF];
goal_config = [1,num_DOF];

fileID = fopen('test_positions.txt','w');
fprintf(fileID,'%6s %12s\n','Start Position                                          ','Goal_Position');
for i = 1 : num_samples_to_generate
    rng(i);
    for i = 1 : num_DOF
      start_config(i) = min_val + rand*(max_val-min_val);
    end
    start_config
    fprintf(fileID,'%6.4f\t',start_config);  
    fprintf(fileID,'\t\t'); 
    for i = 1 : num_DOF
      goal_config(i) = min_val + rand*(max_val-min_val);
    end
    goal_config

    fprintf(fileID,'%6.4f\t',goal_config);   
    fprintf(fileID,'\n'); 
end
fclose(fileID);