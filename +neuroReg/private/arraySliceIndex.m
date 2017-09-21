function index = arraySliceIndex(n,I,d)
index=repmat({':'},1,n);
index{d}=I;
end