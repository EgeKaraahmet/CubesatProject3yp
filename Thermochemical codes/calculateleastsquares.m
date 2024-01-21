for i=1:length(AM)
square(i) = AM(i) * AM(i);
end
Averagearea = sum(square);
Averagearea = Averagearea / length(AM);
Averagearea = sqrt(Averagearea);
