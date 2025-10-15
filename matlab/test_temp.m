clc;clear;

order = [2, 4, 4];
fc = [0.2 0.2 0.3];
type = ["low", "low", "high"];

for i = 1:3
    [b,a] = butter(order(i), fc(i), type(i));
    disp(a)
    disp(b)
    disp("------")
end