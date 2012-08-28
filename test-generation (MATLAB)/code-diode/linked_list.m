
function listObject = linked_list(values)

  data = reshape(values,1,[]);
  listObject = struct('display',@display_list,...
                      'addAfter',@add_element,...
                      'delete',@delete_element);

  function display_list
    %# Displays the data in the list
    disp(data);
  end

  function add_element(values,index)
    %# Adds a set of data values after an index in the list, or at the end
    %#   of the list if the index is larger than the number of list elements
    index = min(index,numel(data));
    data = [data(1:index) reshape(values,1,[]) data(index+1:end)];
  end

  function delete_element(index)
    %# Deletes an element at an index in the list
    data(index) = [];
  end

end