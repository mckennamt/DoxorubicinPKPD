function [drugAIF, Y_bound, time_vector, Y_free] = generateChemoAIF(ExperimentTime, numSamps, ...
    compartmentParms, DrugSchedule)


%keyboard
time_vector = linspace(0,ExperimentTime,numSamps+1);%in hours

drugAIF = zeros(size(time_vector));

drugChanges = DrugSchedule.drugChangeTimes;
drugConc = DrugSchedule.drugConc;

for dAdds = 1:length(drugChanges)
    tp = find(time_vector-drugChanges(dAdds) > 0,1);
    drugAIF(tp:end) = drugConc(dAdds);
end

%%ode solve
%tp = find(time_vector-drugChanges(2) > 0,1);
%S0 = [0 0];
%[td_1, Yd1] = ode45(@(t,y) chemoCompartmentModel_drugOn(t,y,compartmentParms,drugAIF,time_vector),time_vector(1:tp-1),S0);
%S0 = [0 Yd1(end,:)];
%[td_2, Yd2] = ode45(@(t,y) chemoCompartmentModel_drugOff(t,y,compartmentParms),time_vector(tp:end),S0);
%drugAIF(tp:end) = drugAIF(tp:end) + Yd2(:,1)';
%Y_bound = [Yd1(:,2);Yd2(:,3)]';
%Y_free = [Yd1(:,1);Yd2(:,2)]';

%analytic solution, 3 compartment
p = compartmentParms;
p(2) = p(1)/p(2);
Y_bound = modelCt3C3P(p,drugAIF,time_vector);
Y_free = Y_bound;