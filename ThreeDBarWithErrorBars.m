function ThreeDBarWithErrorBars(PeakVals, ErrorBarSizes,B1,minA,maxA)

%{
%to test the function run the following code at the matlab prompt:
cols = 5;       %columns in the matrix
rows = 4;       %rows in the matrix
PeakVals = rand(rows,cols);  % matrix of peak heights
Peak_Stdevs = ones(rows,cols)*0.15; % matrix of standard deviations
ThreeDBarWithErrorBars(PeakVals, Peak_Stdevs)
%}

if size(PeakVals) == size(ErrorBarSizes)

    [rows cols] = size(PeakVals);
    % Plot the 3D bar graph
    c=mesh((PeakVals)) ;
    for i=1:length(c)
        c(i).LineWidth=1;
        c(i).EdgeColor='red';
    end
    colormap('white');
    zlim([minA maxA])
    ylim([1 5]);
    xlim([1 7]);
    %Now plot the standard deviations on top of the bars
    hold on;
%     for row = 1:rows,
%         for col = 1:cols,
%             z = PeakVals(row,col) : (ErrorBarSizes(row,col) / 50) : (PeakVals(row,col) + ErrorBarSizes(row,col));
%             x(1:length(z)) = row;
%             y(1:length(z)) = col;
%             if rows >= cols
%                 plot3(y, x, z,'r-')
%             end
%             if rows < cols
%                 plot3(x, y, z,'r-','linewidth', 1)
%             end
%             clear x; clear y; clear z;
%         end
%     end
%     hold off
%     %disp('done');
    for i = 1:size(PeakVals, 1)

% Traverse the x values, or columns
        for j = 1:size(PeakVals, 2)

            X = [j, j];

            % Set the y coordinates for the line
            Y = [i, i];

            % The endpoint of the line
            z_end2 = PeakVals(i, j) + ErrorBarSizes(i, j);
            z_end1 = PeakVals(i, j) - ErrorBarSizes(i, j);


            % Set the z coordinates for the line
            Z = [PeakVals(i,j), z_end2];
      

            % Draw a solid black line according to its coordinates
            plot3(X, Y, Z, 'k-','linewidth',1);
            X = [j - 0.05, j + 0.05];

% Finally, the Z coordinates (Y does not change)
            Z = [z_end2, z_end2];

% Plot the end marker
            plot3(X,Y,Z,'k-','linewidth',1);
            zlim([minA maxA])
            ylim([1 5]);
            xlim([1 7]);
        end % end if
    end
else
    disp('error, matrices not the same size');
end


