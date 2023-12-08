def largest_common_substring(string1, string2, numbering1, numbering2):
    len1 = len(string1)
    len2 = len(string2)

    # Initialize a 2D array to store lengths of common suffixes
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]

    # Variables to store length of the longest common substring and its ending position
    max_len = 0
    end_pos1 = 0
    end_pos2 = 0

    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            if string1[i - 1] == string2[j - 1] or string2[j - 1] == '_':
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > max_len:
                    max_len = dp[i][j]
                    end_pos1 = i
                    end_pos2 = j
            else:
                dp[i][j] = 0

    # Extract the largest common substring from string2
    common_substring = string2[end_pos2 - max_len : end_pos2]

    # Extract the numbering of the common substring
    numbering_common_substring = numbering2[end_pos2 - max_len : end_pos2]

    # Extract the numbering of the starting position in both strings
    start_numbering1 = numbering1[end_pos1 - max_len]
    start_numbering2 = numbering2[end_pos2 - max_len]

    return (
        common_substring,
        numbering_common_substring,
        start_numbering1,
        start_numbering2
    )

# Test with your example
string1 = "CCGCCGCGCCAUGCCUGUGGCG"
string2 = "CCGCCGCGCCAUGCCUGUGGCG"
numbering1 = list(range(1, len(string1) + 1))
numbering2 = list(range(1, len(string2) + 1))

result, numbering_common_substring, start_numbering1, start_numbering2 = largest_common_substring(
    string1, string2, numbering1, numbering2
)
print(string1)
print(string2)
print(result)

print(f"Largest common substring from string2: {result}")
print(f"Numbering of the common substring: {numbering_common_substring}")
print(f"Numbering of the starting position in string1: {start_numbering1}")
print(f"Numbering of the starting position in string2: {start_numbering2}")
