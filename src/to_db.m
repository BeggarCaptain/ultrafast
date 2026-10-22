function B = to_db(A, dr)
%TO_DB  Envelope -> dB image clipped at -DR relative to the global maximum.

A = double(A);
B = 20*log10(max(A, realmin) / max(A(:)));
if nargin > 1 && ~isempty(dr)
    B = max(B, -dr);
end
end
