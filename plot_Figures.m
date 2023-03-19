figure;
x1=transpose(XL_guess(1,:));
y1=transpose(XL_guess(2,:));
z1=transpose(XL_guess(3,:));
x2=transpose(XL_real(1,:));
y2 = transpose(XL_real(2,:));
z2=transpose(XL_real(3,:));
scatter3(x1,y1 ,z1, 'b', 'filled');
scatter3(x2,y2 ,z2, 'r', 'filled');

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

for i=1:length(x1)
  text(x1(i), y1(i), z1(i), num2str(i), 'FontSize', 1, 'FontWeight', 'bold');
  text(x2(i), y2(i), z2(i), num2str(i), 'FontSize', 1, 'FontWeight', 'bold');
end
axis equal;
title('Landmarks');

saveas(gcf, 'my_figure.png');