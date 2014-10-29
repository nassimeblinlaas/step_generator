Data = csvread('./morisawa.csv');  
fprintf('Reading done\n')
%%
%plot(Data2(:, 13 ), 'r' );
%hold on
%plot(Data2(:, 25 ), 'g' );

% colonne 1 Left  foot
% colonne 2 Right foot
% colonne 3 Com x
% colonne 4 Com y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build contact points data
is_contact = zeros( length(Data) , 2);
contact_coord = zeros(length(Data), 3 * 2);
normal_vectors = zeros(length(Data), 3 * 2);    % three scalars for vector foreach two contact points

for i = 1:length(Data)
    if Data(i, 13) == 0
       is_contact(i, 1) = 1;
       contact_coord(i, 1) = Data(i, 11);    % coord x LF
       contact_coord(i, 2) = Data(i, 12);    % coord y LF
       contact_coord(i, 3) = Data(i, 13);    % coord z LF
       normal_vectors(i, 1) = 0;    % x1
       normal_vectors(i, 2) = 0;    % y1
       normal_vectors(i, 3) = 1;    % z1
    end
    if Data(i, 25) == 0
        is_contact(i, 2) = 1;
        contact_coord(i, 4) = Data(i, 23);    % coord x RF
        contact_coord(i, 5) = Data(i, 24);    % coord y RF
        contact_coord(i, 6) = Data(i, 25);    % coord z RF
        normal_vectors(i, 4) = 0;    % x2
        normal_vectors(i, 5) = 0;    % y2
        normal_vectors(i, 6) = 1;    % z2
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
is_contact;
figure(1)
plot(contact_coord(:, 1), 'r');
hold on
plot(contact_coord(:, 2), 'g');
plot(contact_coord(:, 3), 'b');
hold off
figure(2)
plot(contact_coord(:, 4), 'd');
hold on
plot(contact_coord(:, 5), 'k');
plot(contact_coord(:, 6), 'c');
hold off