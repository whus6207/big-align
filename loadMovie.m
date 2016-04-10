function movie_genre=loadMovie(filename)
    tenMmovies = importdata(filename);
    for i=1:size(tenMmovies)
       movies_data(i,:)= strsplit(tenMmovies{i},'::');
    end
    genres = {'Action','Adventure','Animation','Children''s','Comedy','Crime','Documentary','Drama','Fantasy','Film-Noir','Horror','Musical','Mystery','Romance','Sci-Fi','Thriller','War','Western'};

    movie_genre=zeros(size(movies_data,1),size(genres,1));
    for i=1:size(movies_data,1)
        for j=1:size(genres,2)
           movie_genre(i,j)=(size(strfind(movies_data{i,3},genres{j}),1)>0);
        end
    end
end