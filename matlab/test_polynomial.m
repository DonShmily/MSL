clc;clear;

x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
y = [1.0, 2.0, 0.0, 5.0, 4.0, 3.0];

p1 = polyfit(x, y, 3);
test_x = [0.5, 1.5, 2.5, 3.5];
test_y1 = polyval(p1, test_x);

p2 = polyfit(x, y, 4);
test_x = [0.5, 1.5, 2.5, 3.5];
test_y2 = polyval(p2, test_x);
