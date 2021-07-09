function checkIfFails(handle, msg)

    passed = false;
    try
        handle();
    catch ME %#ok<NASGU>
       passed = true; 
    end

    if passed
        disp...
        (...
            [msg,...
             '[PASSED]'...
            ]...
        );
    else
        error('Test faild!');
    end
end
